#!/usr/bin/env python
# -* coding: utf-8 -*-

__author__ = 'sujiawei'
__email__ = 'su_jiawei163@163.com'
__date__ = '2022-04-01'
__version__ = 'V1.3'

import os
import sys
import json
import gzip
import click
import pandas as pd
from pathlib import Path
from collections import defaultdict
from multiprocessing import Process

Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg


def get_gtid2index(db):
    """ 获取genus taxid对应的bwa index的路径"""
    gtid2index = defaultdict()
    for fna_index in db.glob('*/*fna'):
        gtid = fna_index.name.split('.')[0]
        gtid2index[gtid] = str(fna_index)
    return gtid2index


def bwa_align(query, reference, outfile):
    cmd = f'{cfg.bwa} mem -t2 \'{reference}\' \'{query}\' | samtools view -hF4 | samtools sort - > {outfile}'
    return cmd


@click.command(no_args_is_help=True)
@click.option('-i', '--indir', help='input dir')
@click.option('--db', help='absolute path to bwa database')
@click.option('-o', '--outdir', help='result output dir')
def main(indir, db, outdir):
    """ 对指定文件夹下的序列文件使用bwa mem进行比对 """
    gtid2index = get_gtid2index(Path(db))
    fasta_dir = Path(outdir) / 'Fasta'
    sh_lst_wfp = Path(outdir).joinpath('BWAcmd.list').open('w')

    for fna_gz in fasta_dir.glob('*/*/*fna.gz'):
        is_empty = True
        with gzip.open(fna_gz, 'rb') as rfp:
            for line in rfp:
                if str(line, 'utf8').strip():
                    is_empty = False
                    break
        if is_empty:
            continue
        gtid = fna_gz.name.split('.')[0]
        genus = fna_gz.parent.name
        kingdom = fna_gz.parent.parent.name
        reference = gtid2index.get(gtid)
        bam_outdir = Path(outdir) / 'BwaAlign' / kingdom
        if not bam_outdir.exists():
            os.system(f'mkdir -p {bam_outdir}')
        outbam = bam_outdir / f'{gtid}.bam'
        cmd = bwa_align(str(fna_gz), reference, outbam)
        forth_col = '/'.join(['Fasta', kingdom, genus, fna_gz.name])
        sh_lst_wfp.write(f'{kingdom}\t{genus}\t{indir}\t{forth_col}\t{cmd}\n')
    sh_lst_wfp.close()


if __name__ == '__main__':
    main()
