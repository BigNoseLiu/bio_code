#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'sujiawei'
__email__ = 'su_jiawei163@163.com'
__date__ = '2021-12-21'
__version__ = 'V1.3'

import os
import sys
import json
import click
import warnings
import pandas as pd
from pathlib import Path
from collections import defaultdict, OrderedDict

warnings.filterwarnings('ignore')

###################################################################################
"""                             inputs data types                               """
input_xls_dtypes = OrderedDict((
    ('tax_id',                 int),
    ('species',                str),
    ('reads_count',            int),
    ('coverage',               int),
    ('species_cn',             str),
    ('genus_cn',               str),
    ('noun',                   str),
    ('medicine',               str),
    ('ncbi',                   str),
    ('blast_maxscore',         str),
    ('mark',                   str),
    ('type_judge',             str)
))

"""                             output data types                               """
# 质控结果
sheet_data = OrderedDict((
    ('file',              str),
    ('Num reads',         int),
    ('Num Bases',         int),
))
# 耐药
sheet_ARG = OrderedDict((
    ('#ARG',             str),
    ('Coverage',         int),
    ('antibiotic',       str),
    ('Pathogen',         str),
))
# 毒力
sheet_VF = OrderedDict((
    ('refname',         str),
    ('Gene',            str),
    ('VFgene',          str),
    ('Summary',         str),
    ('Species',         str),
    ('Identity',        float),
))
# 物种鉴定
sheet_Pathogen = OrderedDict((
    ('Kingdom',        str),
    ('GenusName',      str),
    ('GenusCN',        str),
    ('GenusReads',     int),
    ('GroupName',      str),
    ('GroupReads',     int),
    ('ScientificName', str),
    ('ChineseName',    str),
    ('Comment',        str),
    ('Reads_Number',   int),
    ('GMRN',           int),
    ('TPM',            float),
    ('RH-RPM',         float),
    ('Normalization',  int),
    ('Verification',   str),
    ('Focus',          str),
    ('Coverage',       float),
))
###################################################################################


def get_qc_df(in_dir, sample, library):
    # 获取质控前后数据的统计结果
    fastp_smry = Path(in_dir) / "01.QC/01.Fastp" / sample / library / f"{sample}.json"
    with open(fastp_smry, 'r') as fs_fp:
        fastp_stats = json.load(fs_fp)
    data = {
            'file': [f"{sample}_raw", f"{sample}_clean"],
            'Num reads': [fastp_stats['summary']['before_filtering']['total_reads'],
                          fastp_stats['summary']['after_filtering']['total_reads']],
            'Num Bases': [fastp_stats['summary']['before_filtering']['total_bases'],
                          fastp_stats['summary']['after_filtering']['total_bases']]
            }
    try:
        df = pd.DataFrame(data)
    except:
        print(data)
        sys.exit(1)
    df = df.astype(sheet_data)
    return df, data


def get_arg_df(in_dir, sample, library):
    """ 耐压基因检测结果格式转换 """
    in_file = Path(in_dir) / "02.Annotation/02.ResistanceGene" / sample / library / "ARG.cov.pick.xls"
    if os.path.exists(in_file):
        in_df = pd.read_csv(in_file, sep='\t', header=0)
        in_df = in_df.rename(columns={'coverage': 'Coverage', 'pathogen': 'Pathogen'})
        out_df = in_df[sheet_ARG.keys()]
    else:
        cols = list(sheet_ARG.keys())
        out_df = pd.DataFrame(columns=cols)
    return out_df


def get_vf_df(in_dir, sample, library):
    """ 毒力因子检测结果格式转换 """
    in_file = Path(in_dir) / "02.Annotation/03.VFDB" / sample / library / "VF.xls"
    if in_file.exists():
        in_df = pd.read_csv(in_file, sep='\t', header=0)
        in_df = in_df.rename(columns={'gene': 'Gene', 'summary': 'Summary', 'species': 'Species', 'identity': 'Identity'})
        out_df = in_df[sheet_VF.keys()]
    else:
        cols = ['tax_id', 'refname', 'gene', 'VFgene', 'summary', 'species', 'identity', 'type_judge']
        out_df = pd.DataFrame(columns=cols)
    return out_df


def get_tax_df(in_dir, sample):
    """ 病原鉴定结果格式转换 """
    tax_xlsx = Path(in_dir) / f'plugin_out/Result/{sample}/Total_Detail.xlsx'
    df = pd.read_excel(tax_xlsx)
    return df


def generate_summary(qc_stats, tax_df):
    """ 根据fastp质控数据、物种鉴定结果生成summary文件 """
    summary = defaultdict(list)
    summary['NumReads'].append(qc_stats['Num reads'][0])
    summary['NumBases'].append(qc_stats['Num Bases'][0])
    summary['NumSpcReads'].append(sum(tax_df['Reads_Number'].tolist()))
    df_summary = pd.DataFrame(summary)
    return df_summary


@click.command(no_args_is_help=True)
@click.option('-w', '--workdir', help='Input dir, where the PIP outputs are.')
@click.option('-n', '--name', help='The sample name')
@click.option('-l', '--library', help='The library of the sample, DNA or RNA.')
def workflow(workdir, name, library):
    """ 将流程的分析结果按指定格式汇总为一个xls文件 """
    outdir = Path(workdir) / 'plugin_out/Result' / name
    if not Path(outdir).exists():
        os.system(f"mkdir -p {outdir}")

    # merge sheet to xls
    sheet_name2df = OrderedDict()
    sheet_name2df['Data'], qc_stats = get_qc_df(workdir, name, library)
    sheet_name2df['ARG'] = get_arg_df(workdir, name, library)
    sheet_name2df['VF'] = get_vf_df(workdir, name, library)
    sheet_name2df['Pathogen'] = get_tax_df(workdir, name)

    out_xlsx = Path(outdir) / f"{name}.merge.xls"
    with pd.ExcelWriter(out_xlsx) as writer:
        for sheet_name, df in sheet_name2df.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    df_summary = generate_summary(qc_stats, sheet_name2df['Pathogen'])
    out_summary = Path(outdir) / 'raw_summary.xls'
    df_summary.to_csv(out_summary, sep='\t', index=False)


if __name__ == '__main__':
    workflow()
