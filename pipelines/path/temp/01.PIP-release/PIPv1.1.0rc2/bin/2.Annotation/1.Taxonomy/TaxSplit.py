#!/usr/bin/env python
import os
import sys
import gzip
import click
import warnings
import pandas as pd
from pathlib import Path
from multiprocessing import Process, Pool
Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg
from configure import createPath, checkFile, checkType, printInfo, Resource

warnings.filterwarnings('ignore')

############################ developer info ##############################
__AUTHOR__ = "Sujiawei"
__EMAIL__ = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__ = "2022-05-11"
##########################################################################

KINGDOM2TAG = {
    'Bacteria': 'bacteria',
    'Eukaryota:Fungi': 'fungi',
    'Eukaryota:Parasite': 'parasite',
    'Eukaryota:Protozoa': 'protozoa',
    'SpecialPathogens': 'special_pathogen',
    'Viruses': 'virus'
}


def check_file(file):
    if not Path(file).exists():
        raise FileExistsError(f'{str(file)} not exists!')
    return


############################### Annotating by DB ###############################
def get_tax_info(db):
    """ get seq_taxid annotation from PIDB """
    infile1 = Path(db) / 'all.Taxonomy.txt'
    infile2 = Path(db) / 'all.Taxonomy.other.txt'
    infile3 = Path(db) / 'GroupList.info'
    infile4 = Path(db) / 'seqTaxid2S1LevelTaxid'
    for file in infile1, infile2, infile3:
        check_file(file)

    tid2spe_tid = {}
    tid2s1_tid = {}
    tid2kingdom = {}
    tid2name = {}
    gns_tid_set = set()
    grp_tid_set = set()
    spe_tid2grp_tid = {}

    for file in infile1, infile2:
        df = pd.read_csv(file, sep='\t', dtype=str)
        for row in df.itertuples():
            tid2spe_tid[row[1]] = row[2]
            tid2spe_tid[row[2]] = row[2]
            tid2kingdom[row[1]] = row[3]
            tid2kingdom[row[2]] = row[3]
            tid2name[row[1]] = row[10]
            tid2name[row[2]] = row[9]
            tid2name[row[14]] = row[8]
            gns_tid_set.add(row[14])
            tid2kingdom[row[14]] = row[3]

    df = pd.read_csv(infile3, sep='\t', dtype=str)
    for row in df.itertuples():
        grp_tid_set.add(row[1])
        spe_tid2grp_tid[row[3]] = row[1]
        tid2kingdom[row[1]] = row[9]
        tid2name[row[1]] = row[2]

    df = pd.read_csv(infile4, sep='\t', dtype=str)
    for row in df.itertuples():
        tid2s1_tid[row[1]] = row[2]

    return {
        'tid2spe_tid': tid2spe_tid,
        'tid2s1_tid': tid2s1_tid,
        'tid2kingdom': tid2kingdom,
        'tid2name': tid2name,
        'gns_tid_set': gns_tid_set,
        'grp_tid_set': grp_tid_set,
        'spe_tid2grp_tid': spe_tid2grp_tid
    }


def build(in_fq, seq_id2tax_id, outdir, thread):
    """ 根据物种鉴定所得的序列构建覆盖图用的临时数据库 """
    db_path = Path(outdir) / 'Coverage' / '01.db_path'
    dir_lib = db_path / 'library'
    dir_tax = db_path / 'taxonomy'
    for path in dir_lib, dir_tax:
        if not path.exists():
            os.system(f'mkdir -p {str(path)}')
    # output seq
    out_seq = dir_lib / f'library.fna'
    prelim_map = dir_lib / 'prelim_map.txt'

    with gzip.open(in_fq, 'rb') as rfp, open(out_seq, 'w') as wfp_seq, open(prelim_map, 'w') as wfp_map:
        tag = False
        n = 0
        for line in rfp:
            n += 1
            if n % 4 == 1:
                seq_id = str(line, 'utf-8').strip().replace('@', '', 1)
                tax_id = seq_id2tax_id.get(seq_id)
                wfp_seq.write(f'>kraken:taxid|{tax_id}|{seq_id}\n')
                wfp_map.write(f'TAXID\tkraken:taxid|{tax_id}|{seq_id}\t{tax_id}\n')
                tag = True
            elif tag:
                wfp_seq.write(str(line, 'utf-8'))
                tag = False

    # link taxonomy files
    os.system(f'ln -s {cfg.taxonomy}/*dmp {str(dir_tax)}')
    os.system(f'ln -s {cfg.taxonomy}/*taxid {str(dir_tax)}')

    # build_db: OMP only wants you to use 2 threads
    cmd = f'{cfg.kraken2build} --build --db {str(db_path)} --threads 2'
    os.system(cmd)


def seq2df(fq_list, seq_id2tax_id, unique_seqid_set, tax_info):
    """ 将序列文件fq转为df """
    # deal with seq file
    data_seq = []
    tag = False
    for fq_file in fq_list.split(','):
        with gzip.open(fq_file, 'rb') as rfp:
            n = 0
            for line in rfp:
                n += 1
                if n % 4 == 1:
                    ll = str(line, 'utf-8').strip().split(' ')
                    item = [ll[0].replace('@', '', 1)]
                    tag = True
                elif tag:
                    item.append(str(line, 'utf-8').strip())
                    data_seq.append(item)
                    del item
                    tag = False
    df = pd.DataFrame(data_seq, columns=['SeqID', 'Seq'])
    df['Taxid'] = [seq_id2tax_id[i] for i in df['SeqID']]

    # unique mapping
    df_unique = df[df['SeqID'].isin(unique_seqid_set)]

    # unique mapping, species
    df_spe = df_unique[df_unique['Taxid'].isin(tax_info['tid2spe_tid'])]
    df_spe['SpeTaxid'] = [tax_info['tid2spe_tid'][i] for i in df_spe['Taxid']]
    df_spe['Kingdom'] = [tax_info['tid2kingdom'][i] for i in df_spe['Taxid']]

    # unique mapping, group/complex
    df_grp = df[df['Taxid'].isin(tax_info['grp_tid_set'])]
    df_grp['SpeTaxid'] = df_grp['Taxid']
    df_grp['Kingdom'] = [tax_info['tid2kingdom'][i] for i in df_grp['Taxid']]
    df_merge = pd.concat([df_spe, df_grp])

    return df_merge


def split_seq(df_seq, df_res, tax_info, kingdom, out_dir):
    """ 按kingdom分类信息输出unique/multi mapping reads """
    df_kd_res = df_res[df_res['Kingdom'] == kingdom]
    if df_kd_res.empty:
        return
    target_spe = set(df_kd_res['SpeciesSN'])

    df_kd_seq = df_seq[df_seq['Kingdom'] == kingdom]
    if df_kd_seq.empty:
        return
    df_kd_seq = df_kd_seq.sort_values('SpeTaxid')

    seq_file = out_dir / f'{KINGDOM2TAG.get(kingdom)}.result.fa.gz'
    with gzip.open(seq_file, 'w') as wfp_seq:
        if kingdom == 'Viruses':
            for row in df_kd_seq.itertuples():
                if row[3] in tax_info['tid2s1_tid']:
                    spe_name = tax_info['tid2name'].get(tax_info['tid2s1_tid'][row[3]])
                else:
                    spe_name = tax_info['tid2name'].get(row[4])
                if spe_name in target_spe:
                    line = f'>kraken:taxid|{row[3]}|{row[1]} {spe_name}\n{row[2]}\n'
                    wfp_seq.write(bytes(line, 'utf-8'))
        else:
            for row in df_kd_seq.itertuples():
                spe_name = tax_info['tid2name'].get(row[4])
                if spe_name in target_spe:
                    line = f'>kraken:taxid|{row[3]}|{row[1]} {spe_name}\n{row[2]}\n'
                    wfp_seq.write(bytes(line, 'utf-8'))


@click.command(no_args_is_help=True)
@click.option('-k', '--kraken_out', help='Kraken2 output file <str>')
@click.option('-r', '--result', help='Pathogens taxonomy identification result, xlsx. <str>')
@click.option('--db', help='Path of PIDB <str>')
@click.option('--fq', help='Clean data\'s fastq.gz file, separated with comma, example: --fq fq1,fq2 <str>')
@click.option('-o', '--outdir', help='Results output dir <str>')
@click.option('-t', '--thread', help='Specify how many threads can be used.')
def main(kraken_out, result, db, fq, outdir, thread):
    """ Split sequence by tax """
    # check file
    for file in kraken_out, result:
        check_file(file)
    for file in fq.split(','):
        check_file(file)

    # get tax id annotation from kraken output
    seq_id2tax_id = {}
    unique_seqid_set = set()
    with open(kraken_out, 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.split('\t')
            seq_id2tax_id[ll[1]] = ll[2]
            unique_seqid_set.add(ll[1])

    process_lst = []

    # build kraken db for classified seq
    p = Process(target=build, args=(fq, seq_id2tax_id, Path(outdir), thread))
    p.start()
    process_lst.append(p)

    # get tax info from PIDB
    tax_info = get_tax_info(db)
    # get df of sequences with info
    df_seq = seq2df(fq, seq_id2tax_id, unique_seqid_set, tax_info)
    if df_seq.empty:
        print(f'{fq} is empty. Please check it.')
        return

    # get df of taxonomy result
    df_res = pd.read_excel(result)
    if df_res.empty:
        print(f'{result} is empty. Please check it.')
        return

    outdir = Path(outdir) / 'ResultSeq'
    if not outdir.exists():
        os.system(f'mkdir -p {str(outdir)}')

    # split seq
    for kd in KINGDOM2TAG.keys():
        p = Process(target=split_seq, args=(df_seq, df_res, tax_info, kd, outdir))
        p.start()
        process_lst.append(p)

    for p in process_lst:
        p.join()


if __name__ == '__main__':
    main()