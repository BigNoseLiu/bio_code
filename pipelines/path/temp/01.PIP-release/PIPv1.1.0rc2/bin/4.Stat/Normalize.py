#!/usr/bin/env python
import os
import sys
import json
import math
import click
import pandas as pd
from pathlib import Path

############################ developer info ##############################
__AUTHOR__ = "Sujiawei"
__EMAIL__ = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__ = "2022-05-06"
##########################################################################


@click.command(no_args_is_help=True)
@click.option('-w', '--workdir', help='Tasks working directory')
@click.option('-s', '--samples', help='Sample list')
@click.option('-n', '--n_million', help='N million to normalize [20]', default=20)
def main(workdir, samples, n_million):
    """ 对检出reads进行n M均一化 """
    qc_dir = Path(workdir) / '01.QC/01.Fastp'
    tax_dir = Path(workdir) / '02.Annotation/01.Taxonomy'
    plg_dir = Path(workdir) / 'plugin_out/Result'

    with open(samples, 'r') as rfp:
        smp2lib = {line.split('\t')[0]: line.split('\t')[1] for line in rfp.readlines()}
    for smp, lib in smp2lib.items():
        qc_jsn = qc_dir / smp / lib / f'{smp}.json'
        with open(qc_jsn, 'r') as rfp:
            qc_info = json.load(rfp)
        coe = n_million * 10 ** 6 / qc_info["summary"]["before_filtering"]["total_reads"]
        total_file = tax_dir / smp / lib / 'Total_Detail.xlsx'
        patho_file = tax_dir / smp / lib / 'Pathogeny_Detail.xlsx'
        out_dir = plg_dir / smp
        if not out_dir.exists():
            os.system(f'mkdir -p {str(out_dir)}')
        for file in total_file, patho_file:
            df = pd.read_excel(file)
            df['GenusReads(Normalization)'] = df.apply(lambda x: math.ceil(x['GenusReads'] * coe), axis=1)
            df['GroupReads(Normalization)'] = df.apply(lambda x: math.ceil(x['GroupReads'] * coe), axis=1)
            df['SpeciesReads(Normalization)'] = df.apply(lambda x: math.ceil(x['SpeciesReads'] * coe), axis=1)
            df = df.sort_values(by='SpeciesReads', ascending=False)
            df_out = pd.DataFrame(columns=df.columns)
            # sorted by kingdom, genus reads, species reads
            for kd in ['Bacteria', 'SpecialPathogens', 'Eukaryota:Fungi', 'Eukaryota:Parasite',
                       'Eukaryota:Protozoa', 'Viruses']:
                df_kd = df[df['Kingdom'] == kd]
                if not df_kd.empty:
                    df_kd = df_kd.sort_values(by='SpeciesReads', ascending=False)
                    gns_dct = {gns: 0 for gns in df_kd['GenusSN']}
                    for gns in gns_dct.keys():
                        df_temp = df_kd[df_kd['GenusSN'] == gns]
                        df_temp = df_temp.sort_values(by='SpeciesReads', ascending=False)
                        df_out = pd.concat([df_out, df_temp])

            # output
            df_out.to_csv(out_dir / file.with_suffix('.txt').name, index=False, sep='\t')
            df_out.to_excel(out_dir / file.name, index=False)


if __name__ == '__main__':
    main()