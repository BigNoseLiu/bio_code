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
__DATE__ = "2022-06-01"
##########################################################################

KINGDOM2TAG = {
    'Bacteria': 'Bacteria',
    'Eukaryota:Fungi': 'Eukaryota_Fungi',
    'Eukaryota:Parasite': 'Eukaryota_Parasite',
    'Eukaryota:Protozoa': 'Eukaryota_Protozoa',
    'SpecialPathogens': 'SpecialPathogens',
    'Viruses': 'Viruses'
}


def check_file(file):
    if not Path(file).exists():
        raise FileExistsError(f'{str(file)} not exists!')
    return


############################### Annotating by DB ###############################
def get_tax_info(db):
    """ 获取数据库中物种的注释信息 """
    tax_file = Path(db) / 'all.Taxonomy.txt'
    check_file(tax_file)

    tid2name = {}
    tid2spe_tid = {}
    tid2kingdom = {}
    gns_tid_set = set()

    df = pd.read_csv(tax_file, sep='\t', dtype=str)
    for row in df.itertuples():
        tid2spe_tid[row[3]] = row[4]
        tid2kingdom[row[3]] = row[5]
        tid2kingdom[row[4]] = row[5]
        gns_tid_set.add(row[17])
        if row[5] == 'Viruses':
            tid2name[row[3]] = row[13]
            tid2name[row[4]] = row[13]
        else:
            tid2name[row[3]] = row[11]
            tid2name[row[4]] = row[11]

    # Group infos for target groups
    grp_file = Path(db) / 'GroupList.info'
    check_file(grp_file)
    df_grp = pd.read_csv(grp_file, sep='\t', dtype=str)
    df_grp = df_grp[df_grp['Focus'] == '*']
    for row in df_grp.itertuples():
        tid2name[row[1]] = row[2]
        tid2kingdom[row[1]] = row[9]

    grp_tid_set = set(df_grp['GroupTaxid'])
    spe_tid2grp_tid = {row[1]: row[2] for row in df_grp[['SpeciesTaxid', 'GroupTaxid']].itertuples()}

    return {
        'tid2name': tid2name,
        'tid2spe_tid': tid2spe_tid,
        'tid2kingdom': tid2kingdom,
        'gns_tid_set': gns_tid_set,
        'grp_tid_set': grp_tid_set,
        'spe_tid2grp_tid': spe_tid2grp_tid,
    }


def build(in_fq, seq_id2tax_id, outdir, db, thread):
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
    taxonomy_dir = Path(db) / 'taxonomy'
    os.system(f'ln -s {str(taxonomy_dir)}/*dmp {str(dir_tax)}')
    os.system(f'ln -s {str(taxonomy_dir)}/*taxid {str(dir_tax)}')

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

    seq_file = out_dir / f'{KINGDOM2TAG.get(kingdom)}.result.fq.gz'
    with gzip.open(seq_file, 'w') as wfp_seq:
        for row in df_kd_seq.itertuples():
            spe_name = tax_info['tid2name'].get(row[3])
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

    outdir = Path(outdir)

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
    p = Process(target=build, args=(fq, seq_id2tax_id, outdir, db, thread))
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