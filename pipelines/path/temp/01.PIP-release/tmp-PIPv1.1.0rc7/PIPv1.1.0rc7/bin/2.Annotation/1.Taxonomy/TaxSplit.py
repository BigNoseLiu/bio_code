#!/usr/bin/env python
import os
import sys
import gzip
import time

import click
import warnings
import pandas as pd
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict

Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg
from configure import createPath, checkFile, checkType, printInfo, Resource, LogExceptions

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
            # tid2name[row[4]] = row[13]
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
    spe2grp = {row[1]: row[2] for row in df_grp[['SpeciesName', 'GroupName']].itertuples()}

    return {
        'tid2name': tid2name,
        'tid2spe_tid': tid2spe_tid,
        'tid2kingdom': tid2kingdom,
        'gns_tid_set': gns_tid_set,
        'grp_tid_set': grp_tid_set,
        'spe_tid2grp_tid': spe_tid2grp_tid,
        'spe2grp': spe2grp,
    }


def get_seqid2seq(infile):
    seq_id2seq = {}
    if str(infile).endswith('gz'):
        with gzip.open(infile, 'rb') as rfp:
            n = 0
            tag = False
            seq_id = ''
            for line in rfp:
                n += 1
                if n % 4 == 1:
                    seq_id = str(line, 'utf-8').strip().split(' ')[0].replace('@', '', 1)
                    tag = True
                elif tag:
                    seq_id2seq[seq_id] = str(line, 'utf-8')
                    tag = False
                    seq_id = ''
    else:
        with open(infile, 'r') as rfp:
            tag = False
            seq_id = ''
            for n, line in enumerate(rfp, 1):
                if n % 4 == 1:
                    seq_id = line.strip().split(' ')[0].replace('@', '', 1)
                    tag = True
                elif tag:
                    seq_id2seq[seq_id] = line
                    tag = False
                    seq_id = ''
    return seq_id2seq


def build(seq_id2seq, seq_id2tax_id, outdir, db, thread):
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

    with open(out_seq, 'w') as wfp_seq, open(prelim_map, 'w') as wfp_map:
        for seq_id, seq in seq_id2seq.items():
            tax_id = seq_id2tax_id.get(seq_id, '0')
            wfp_seq.write(f'>kraken:taxid|{tax_id}|{seq_id}\n{seq}')
            wfp_map.write(f'TAXID\tkraken:taxid|{tax_id}|{seq_id}\t{tax_id}\n')

    # link taxonomy files
    taxonomy_dir = Path(db) / 'taxonomy'
    os.system(f'ln -s {str(taxonomy_dir)}/*dmp {str(dir_tax)}')
    os.system(f'ln -s {str(taxonomy_dir)}/*taxid {str(dir_tax)}')

    # build_db: OMP only wants you to use 2 threads
    cmd = f'{cfg.kraken2build} --build --db {str(db_path)} --threads 2'
    os.system(cmd)


def split_seq(seq_id2seq, seq_id2tax_id, target_spe2kd, tax_info):
    """ 输出物种序列 """
    spe2seq = defaultdict(list)

    for seq_id, seq in seq_id2seq.items():
        tax_id = seq_id2tax_id.get(seq_id)
        if not tax_id:
            continue
        spe_name = tax_info['tid2name'].get(tax_id)
        grp_name = tax_info['spe2grp'].get(spe_name)
        if spe_name in target_spe2kd :
            if grp_name:
                # 复合群
                line = ''.join(['>', seq_id, ' ', grp_name, ' | ', spe_name, '\n', seq])
            else:
                line = ''.join(['>', seq_id, ' ', spe_name, '\n', seq])
            spe2seq[spe_name].append(bytes(line, 'utf-8'))

    return spe2seq


def write_data(spe2data, outfile):
    with gzip.open(outfile, 'wb') as wfp:
        for seq_lst in spe2data.values():
            wfp.write(b''.join(seq_lst))


@click.command(no_args_is_help=True)
@click.option('-k', '--kraken_out', help='Kraken2 output file <str>')
@click.option('-r', '--result', help='Pathogens taxonomy identification result, xlsx. <str>')
@click.option('--db', help='Path of PIDB <str>')
@click.option('--fq', help='Clean data\'s fastq.gz file, separated with comma, example: --fq fq1,fq2 <str>')
@click.option('-o', '--outdir', help='Results output dir <str>')
@click.option('-t', '--thread', help='Specify how many threads can be used.')
def main(kraken_out, result, db, fq, outdir, thread):
    """ Split sequence by tax """
    start_time = time.time()

    # check file
    for file in kraken_out, result:
        check_file(file)
    for file in fq.split(','):
        check_file(file)

    outdir = Path(outdir)
    if not outdir.exists():
        os.system(f'mkdir -p {str(outdir)}')

    # get tax id annotation from kraken output
    seq_id2tax_id = {}
    with open(kraken_out, 'r') as rfp:
        for _, line in enumerate(rfp):
            if line.startswith('>') or '\t9606\t' in line:
                continue
            ll = line.split('\t')
            seq_id2tax_id[ll[1]] = ll[2]

    # get seq from input
    seq_id2seq = get_seqid2seq(fq)
    if not seq_id2seq:
        print(f'{fq} is empty. Please check it.')
        return

    if int(thread) > 2:
        parallel_num = int(thread) - 2
    else:
        parallel_num = int(thread)
    pool = Pool(parallel_num)

    # build kraken db for classified seq
    pool.apply_async(LogExceptions(build), args=(seq_id2seq, seq_id2tax_id, outdir, db, thread))

    # get df of taxonomy result
    df_res = pd.read_csv(result, sep='\t', dtype=str)
    if df_res.empty:
        print(f'{result} is empty. Please check it.')
        return
    target_spe2kd = {}
    for row in df_res[['Kingdom', 'SpeciesSN']].itertuples():
        target_spe2kd[row[2]] = row[1]

    # get tax info from PIDB
    tax_info = get_tax_info(db)

    # split seq
    spe2seq_lst = split_seq(seq_id2seq, seq_id2tax_id, target_spe2kd, tax_info)

    # sort by species
    kd2data = {kd: defaultdict(list) for kd in set(df_res['Kingdom'])}
    for spe, spe2seq_lst in spe2seq_lst.items():
        kd = target_spe2kd.get(spe)
        try:
            kd2data[kd][spe].extend(spe2seq_lst)
        except:
            print(kd, spe)
            sys.exit(1)

    for kd in kd2data.keys():
        outfile = outdir / f'{KINGDOM2TAG.get(kd)}.result.fq.gz'
        pool.apply_async(LogExceptions(write_data), args=(kd2data.get(kd), outfile))
    pool.close()
    pool.join()

    print(f'TaxSplit Elapsed Time: ', time.time() - start_time)


if __name__ == '__main__':
    main()
