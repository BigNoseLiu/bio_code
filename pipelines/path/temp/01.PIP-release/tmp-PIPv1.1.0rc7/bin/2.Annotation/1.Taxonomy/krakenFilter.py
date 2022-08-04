#!/usr/bin/env python
# -* coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2022-07-01'
__version__ = 'V1.0'

import os
import sys
import fileinput
from sys import argv
from argparse import ArgumentParser

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../conf'))
import configure as cfg
from configure import checkFile

import click
import pandas as pd
from pathlib import Path


@click.command(no_args_is_help=True)
@click.option('--ko', help="Input kraken's output. <str>")
@click.option('--db', help="Database version of PIP. <str>")
@click.option('--fc', help="Filter criteria. {min, highconf, unique}", default='unique')
@click.option('--out', help="Output file. <str>")
def main(ko, db, fc, out):
    """"""
    checkFile(ko)

    db_tax_file = Path(cfg.KrakenDB) / db / 'all.Taxonomy.txt'
    df = pd.read_csv(db_tax_file, sep='\t')
    tid2spe_tid = {}
    for row in df[['Taxid', 'SpeciesTaxid']].itertuples():
        tid2spe_tid[row[1]] = row[2]

    if fc == 'unique':
        with open(ko, 'r') as rfp, open(out, 'w') as wfp:
            for _, line in enumerate(rfp):
                if not line.startswith('C'):
                    wfp.write(line)
                ll = line.split('\t')
                if ll[2] not in tid2spe_tid:
                    wfp.write(line)
                elif '9606:' in ll[4]:
                    wfp.write(f'C\t{ll[1]}\t131567\t{ll[3]}\t{ll[4]}')
                else:
                    spe_tid_set = set()
                    for detail in ll[4].strip().split(' '):
                        try:
                            key, value = detail.split(':')
                        except:
                            print(detail)
                            sys.exit(1)
                        if key != '|' and key in tid2spe_tid:
                            spe_tid_set.add(tid2spe_tid.get(key))
                    if len(spe_tid_set) > 1:
                        wfp.write(f'C\t{ll[1]}\t131567\t{ll[3]}\t{ll[4]}')
                    else:
                        wfp.write(line)


if __name__ == '__main__':
    main()
