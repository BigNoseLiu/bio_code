#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'sujiawei'
__email__ = 'su_jiawei163@163.com'
__date__ = '2022-03-31'
__version__ = 'V1.3'

import os
import sys
import click
import warnings
import pandas as pd
from pathlib import Path
from collections import defaultdict
from multiprocessing import Process

Bin = Path(__file__).resolve().parent
sys.path.append(os.path.join(Bin, '../../conf'))
import configure as cfg
from configure import checkFile, createPath, printInfo

warnings.filterwarnings('ignore')

# scripts
QC            = Bin / 'QC.stat.pl'
Normalize     = Bin / 'Normalize.py'
CommonSpecies = Bin / 'CommonSpecies.pl'
MergeResults  = Bin / 'mergeResults.py'
FillCoverage  = Bin / 'fillCoverage.py'
SplitQC       = Bin / 'split_QC_Stat.py'


def check_dir(path):
    if not Path(path).exists():
        os.system(f'mkdir -p {str(path)}')


def cp(src, tgt):
    os.system(f'cp -fr {str(src)} {str(tgt)} &')


def execute(cmd):
    os.system(cmd)


@click.command(no_args_is_help=True)
@click.option('-s', '--samples', help='samples list, Format: \"sample DNA|RNA read1.fq [read2.fq]\"\n'
                                      'There are 4 columns separated with Tab in the file, the 1st is \n'
                                      'sample name, the 2nd is sample type (DNA or RNA), the 3rd is the \n'
                                      'fastq file of reads (read1 for PE), the 4th is the fastq file of \n'
                                      'read2 (only for PE).')
@click.option('-a', '--analysis', help='The analytic steps you done. The different steps separated by \'-\', \n'
                                       'the arguments for every step follow the \':\' and separated by \',\'. \n'
                                       'The selectable steps as bollow: default:\n'
                                       '[Merge-QC:Fastp,RemoveHost,rRNAFilter-Annotation:Taxonomy,ResistanceGene,VFDB]')
@click.option('-w', '--workdir', help='The path including 00.Rawdata,01.QC,.... [./]')
def main(samples, analysis, workdir):
    workdir = Path(workdir)

    result_dir = workdir / 'plugin_out/Result'
    check_dir(result_dir)
    qc_dir  = workdir / '01.QC'
    tax_dir = workdir / '02.Annotation/01.Taxonomy'
    arg_dir = workdir / '02.Annotation/02.ResistanceGene'
    vf_dir  = workdir / '02.Annotation/03.VFDB'

    db_ver = analysis.split('=')[1].split(',')[0]
    db_path = Path(cfg.KrakenDB) / db_ver

    os.system(f'{cfg.perl} {QC} {str(qc_dir)}/sample.list {str(result_dir)}/QC.Stat')
    df_qc = pd.read_csv(result_dir / 'QC.Stat', sep='\t', dtype=str)
    df_qc.to_excel(result_dir / 'QC.Stat.xlsx', index=False)
    os.system(f'{cfg.python} {SplitQC} -d {str(workdir)}')
    os.system(f'{cfg.python} {Normalize} --qc {str(result_dir)}/QC.Stat --anno {str(tax_dir)}/sample.list')
    os.system(f'{cfg.perl} {CommonSpecies} --list {str(tax_dir)}/sample.list --out {str(result_dir)}')
    os.chdir(str(result_dir))
    os.system(f'zip -r {str(result_dir)}/zhikong.zip Total_CommonSpecies* Pathogeny_CommonSpecies* QC.Stat*')

    df = pd.read_csv(samples, sep='\t', dtype=str, header=None)

    # basecaller_results
    basecall_dir = workdir / 'basecaller_results'
    check_dir(basecall_dir)
    for row in df.itertuples():
        slide_dir = Path(row[3]).parent
        cp(slide_dir / '*.summaryReport.html', basecall_dir)

    # create by2.run
    with open(workdir / 'by2.run', 'w') as wfp:
        wfp.write('\n'.join(['Result'] + df[0].tolist()))

    process_lst = []
    seq_cp_list = []
    kingdom = ['bacteria', 'fungi', 'LY', 'protozoa', 'viral']
    for row in df.itertuples():
        smp_dir = result_dir / row[1]
        check_dir(smp_dir)
        cmd1 = f'{cfg.python} {FillCoverage} --db {str(db_path)} --indir {str(workdir)} --name {row[1]} ' \
               f'--library {row[2]} --outdir {str(smp_dir)}'
        p = Process(target=execute, args=(cmd1,))
        p.start()
        process_lst.append(p)
        # copy
        orig_tax_dir = tax_dir / row[1] / row[2]
        for kd in kingdom:
            p = Process(target=cp, args=(orig_tax_dir / f'{kd}.xls', smp_dir / f'{kd}.xls'))
            p.start()
            process_lst.append(p)
            p = Process(target=cp, args=(orig_tax_dir / f'{kd}.xls', smp_dir / f'{kd}.pick.xls'))
            p.start()
            process_lst.append(p)
        # ARG
        p = Process(target=cp, args=(arg_dir / row[1] / row[2] / 'ARG.cov.pick.xls', smp_dir))
        p.start()
        process_lst.append(p)
        # VF
        p = Process(target=cp, args=(vf_dir / row[1] / row[2] / 'VF.xls', smp_dir))
        p.start()
        process_lst.append(p)
        # merge pathogen results
        cmd2 = f'{cfg.python} {MergeResults} --workdir {str(workdir)} --name {row[1]} --library {row[2]}'
        print(cmd2)
        p = Process(target=execute, args=(cmd2,))
        p.start()
        process_lst.append(p)
        # by2.sample
        p = Process(target=execute, args=(f'touch {str(workdir)}/by2.sample',))
        p.start()
        process_lst.append(p)
        # execute process
        for p in process_lst:
            p.join()
        process_lst.clear()
        # copy sequences
        for kd in kingdom:
            p = Process(target=cp, args=(tax_dir / row[1] / row[2] / f'{kd}.result.fq.gz',
                                         smp_dir / f'{row[1]}.{kd}.result.fq.gz'))
            p.start()
            seq_cp_list.append(p)
        # copy summary
        for file in ['Pathogeny_Detail', 'Total_Detail', 'Taxonomy_Summary']:
            p = Process(target=cp, args=(tax_dir / row[1] / row[2] / f'{file}.*', smp_dir))
            p.start()
            process_lst.append(p)
        for p in process_lst:
            p.join()
        process_lst.clear()

    for p in seq_cp_list:
        p.join()

    os.system(f'touch {str(workdir)}/by2.done')


if __name__ == '__main__':
    main()