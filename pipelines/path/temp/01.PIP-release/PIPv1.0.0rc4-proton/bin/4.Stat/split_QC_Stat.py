#!/usr/bin/env python
import json
import os
import re
import sys
import click
import pandas as pd
from pathlib import Path

############################ developer info ##############################
__AUTHOR__  = "Sujiawei"
__EMAIL__   = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__    = "2022-03-12"
##########################################################################


@click.command(no_args_is_help=True)
@click.option('-d', '--workdir', help='分析结果路径')
def main(workdir):
    infile = Path(workdir) / 'plugin_out/Result/QC.Stat.xlsx'
    df = pd.read_excel(infile, dtype=str, index_col=0).fillna('')
    for sample in df.columns:
        df.loc['Raw Reads', sample] = str(int(int(df.loc['Raw Reads', sample]) / 1000000)) + ' M'
    for sample in df.columns:
        smp_outfile = Path(workdir) / f'plugin_out/Result/{sample}/{sample}.QC.Stat.xlsx'
        df[sample].to_excel(smp_outfile)

    outfile = infile.with_suffix('.tmp.xlsx')
    df.to_excel(outfile)

    os.system(f'mv {str(outfile)} {str(infile)}')


if __name__ == '__main__':
    main()
