#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'sujiawei'
__email__ = 'su_jiawei163@163.com'
__date__ = '2022-05-31'
__version__ = 'V1.1'

import os
import sys
import click
import warnings
import pandas as pd
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

Bin = Path(__file__).resolve().parent
sys.path.append(os.path.join(Bin, '../../conf'))
import configure as cfg
from configure import checkFile, createPath, printInfo

warnings.filterwarnings('ignore')

# scripts
QC            = Bin / 'QC.stat.py'
Normalize     = Bin / 'Normalize.py'
CommonSpecies = Bin / 'CommonSpecies.pl'
MergeResults  = Bin / 'mergeResults.py'
Supplement    = Bin / 'Supplement.py'


def check_dir(path):
    if not Path(path).exists():
        os.system(f'mkdir -p {str(path)}')


def cp(src, tgt):
    """ execute Linux command: cp """
    if src.is_dir() or src.is_file():
        if not src.exists():
            return
    else:
        l = [i for i in src.parent.glob(src.name)]
        if len(l) == 0:
            return
    os.system(f'cp -fr {str(src)} {str(tgt)}')


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
@click.option('-t', '--thread', help='Specified how many threads can be used. [50]', default=50)
def main(samples, analysis, workdir, thread):
    workdir = Path(workdir)
    df = pd.read_csv(samples, sep='\t', dtype=str, header=None)

    result_dir = workdir / 'plugin_out/Result'
    check_dir(result_dir)
    tax_dir = workdir / '02.Annotation/01.Taxonomy'
    arg_dir = workdir / '02.Annotation/02.ResistanceGene'
    vf_dir  = workdir / '02.Annotation/03.VFDB'

    db_ver = analysis.split('=')[1].split(',')[0]
    db_path = Path(cfg.KrakenDB) / db_ver

    # basecaller_results
    basecall_dir = workdir / 'basecaller_results'
    check_dir(basecall_dir)
    for row in df.itertuples():
        slide_dir = Path(row[3]).parent
        cp(slide_dir / '*.summaryReport.html', basecall_dir)

    # create by2.run
    with open(workdir / 'by2.run', 'w') as wfp:
        wfp.write('\n'.join(['Result'] + df[0].tolist()))

    # QC stat
    os.system(f'{cfg.python} {QC} -w {workdir} -s {samples}')
    
    # # Results Normalization
    # os.system(f'{cfg.python} {Normalize} -w {workdir} -s {samples}')

    # os.system(f'{cfg.perl} {CommonSpecies} --list {str(tax_dir)}/sample.list --out {str(result_dir)}')

    pool = Pool(int(thread))

    for row in df.itertuples():
        smp_dir = result_dir / row[1]
        check_dir(smp_dir)

        # copy taxonomy results
        for file in ['Pathogeny_Detail', 'Total_Detail', 'Taxonomy_Summary']:
            cp(tax_dir / row[1] / row[2] / f'{file}.*', smp_dir)
        # copy ARG results
        pool.apply_async(cp, args=(arg_dir / row[1] / row[2] / 'ARG.cov.pick.*', smp_dir))
        # copy VF results
        pool.apply_async(cp, args=(vf_dir / row[1] / row[2] / 'VF.*', smp_dir))
    pool.close()
    pool.join()

    # Err Temp Convert
    for row in df.itertuples():
        smp_dir = result_dir / row[1]
        check_dir(smp_dir)
        arg_txt = smp_dir / 'ARG.cov.pick.txt'
        arg_xls = arg_txt.with_suffix('.xls')
        os.system(f'cp {str(arg_txt)} {str(arg_xls)}')
        vf_txt = smp_dir / 'VF.txt'
        vf_xls = vf_txt.with_suffix('.xls')
        os.system(f'cp {str(vf_txt)} {str(vf_xls)}')

    # Supplement results
    os.system(f'{cfg.python} {Supplement} -w {workdir} -s {samples} --db {str(db_path)} -t {str(thread)}')

    # merge pathogen results
    pool_merge = Pool(int(thread))
    for row in df.itertuples():
        cmd = f'{cfg.python} {MergeResults} --workdir {str(workdir)} --name {row[1]} --library {row[2]}'
        pool_merge.apply_async(os.system, args=(cmd,))
    pool_merge.close()
    pool_merge.join()

    # pack up results
    os.chdir(str(result_dir))
    os.system(f'zip -r {str(result_dir)}/zhikong.zip Total_CommonSpecies* QC.Stat*')

    # touch file for LIMS pulling
    os.system(f'touch {str(Path(workdir))}/by2.sample')

    pool = Pool(int(thread))
    for row in df.itertuples():
        smp_dir = result_dir / row[1]
        # copy sequences
        pool.apply_async(cp, args=(tax_dir / row[1] / row[2] / '*.result.fq.gz', smp_dir))
        # copy coverage plot
        plot_dir = tax_dir / row[1] / row[2] / 'Coverage/04.plot'
        for kd_dir in plot_dir.glob('*'):
            pool.apply_async(cp, args=(kd_dir, result_dir / row[1] / f'{kd_dir.name}_readscount_coverage_plot'))
    pool.close()
    pool.join()

    # Stat done
    os.system(f'touch {str(workdir)}/by2.done')


if __name__ == '__main__':
    main()