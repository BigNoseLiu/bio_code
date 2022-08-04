#!/usr/bin/env python
import os
import re
import sys
import gzip
import json
import click
import pandas as pd
from pathlib import Path
from collections import defaultdict
Bin = Path(__file__).resolve().parent
sys.path.append(os.path.join(Bin, '../conf'))
import configure as cfg
from configure import checkFile, createPath, parseVersion

############################ developer info ##############################
__AUTHOR__  = "Sujiawei"
__EMAIL__   = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__    = "2022-05-11"
##########################################################################


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
    for file in infile1, infile2, infile3:
        check_file(file)

    tid2spe_tid = {}
    tid2kingdom = {}
    tid2name = {}
    gns_tid_set = set()
    grp_tid_set = set()

    for file in infile1, infile2:
        df = pd.read_csv(file, sep='\t', dtype=str)
        for row in df.itertuples():
            tid2spe_tid[row[1]] = row[2]
            tid2spe_tid[row[2]] = row[2]
            tid2kingdom[row[1]] = row[3]
            tid2kingdom[row[2]] = row[3]
            tid2name[row[1]] = row[9]
            tid2name[row[2]] = row[8]
            gns_tid_set.add(row[14])

    df = pd.read_csv(infile3, sep='\t', dtype=str)
    for row in df.itertuples():
        grp_tid_set.add(row[1])

    return {
        'tid2spe_tid': tid2spe_tid,
        'tid2kingdom': tid2kingdom,
        'gns_tid_set': gns_tid_set,
        'grp_tid_set': grp_tid_set
    }


def write_data(df, spe_tid, spe_name, out_dir):
    """ 按spe_tid输出序列 """
    dir_lib = Path(out_dir) / 'library'
    dir_tax = Path(out_dir) / 'taxonomy'
    for path in dir_lib, dir_tax:
        if not path.exists():
            os.system(f'mkdir -p {str(path)}')

    seq_file = out_dir / f'{spe_tid}.fna'
    prelim_map = out_dir / 'prelim_map'

    with open(seq_file, 'w') as wfp_seq, open(prelim_map, 'w') as wfp_pm:
        df_temp = df[df['SpeTaxid'] == spe_tid]
        for row in df_temp.itertuples():
            wfp_seq.write()

    

@click.command(no_args_is_help=True)
# @click.option('-r', '--report', help='Kraken2\'s report file <str>')
# @click.option('-k', '--kraken_out', help='Kraken2 output file <str>')
# @click.option('-m', '--mark', help='Enter a list of sequences id to be validated. <str>')
# @click.option('--db', help='Path of PIDB <str>')
@click.option('--fq', help='Clean data\'s fastq file, separated with comma, example: --fq fq1,fq2 <str>')
@click.option('-o', '--outdir', help='Results output dir <str>')
def main(report, kraken_out, mark, db, fq, outdir):
    """ Split sequence by tax """
    # check file
    for file in report, kraken_out, mark:
        check_file(file)
    for file in fq.split(','):
        check_file(file)

    # deal with seq file
    data_seq = [['SeqID', 'Taxid', 'Seq']]
    with open(fq, 'r') as rfp:
        for n, line in enumerate(rfp, 1):
            if n % 4 == 1:
                ll = line.strip().split(' ')
                seqid = ll[0].replace('@', '', 1)
                taxid = ll[1].split('|')[1]
                item = [seqid, taxid]
            elif n % 4 == 2:
                item.append(line.strip())
                data_seq.append(item)
                del item

    df = pd.DataFrame(data_seq)
    del data_seq

    # get tax info from PIDB
    tax_info = get_tax_info(db)

    dir_unique = Path(outdir) / 'Coverage' / 'UniqueMapping'

    df_spe = df[df['Taxid'].isin(tax_info['tid2spe_tid'])]
    df_spe['SpeTaxid'] = df_spe.apply(lambda x: tax_info['tid2spe_tid'][x['Taxid']], axis=1)
    for spe_tid in df['SpeTaxid']:
        spe_name = tax_info['']
        seq_dir = dir_unique / tax_info['tid2kingdom'][spe_tid] / spe_tid / 'library'



    for row in df_spe.itertuples():
        spe_tid = tax_info['tid2spe_tid'][row[2]]
        seq_dir = dir_unique / tax_info['tid2spe_tid'][row[2]]
        if not seq_outdir.exists():
            os.system(f'mkdir -p {str(seq_outdir)}')
        seq_file = seq_dir / f'{spe_tid}.fna'


    dir_multi = Path(outdir) / 'Coverage' / 'MultiMapping'




if __name__ == '__main__':
    main()