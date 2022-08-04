#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'sujiawei'
__email__ = 'su_jiawei163@163.com'
__date__ = '2021-12-21'
__version__ = 'V1.3'

import os
import sys
import json
import gzip
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
    ('file',             str),
    ('Num reads',        int),
    ('Num Bases',        int),
))
# 耐药
sheet_ARG = OrderedDict((
    ('#ARG',             str),
    ('Identity(%)',      str),
    ('Coverage(%)',      str),
    ('Average_Depth',    str),
    ('Reads_Num',        str),
    ('Drug_class',       str),
    ('AMR_Family_class', str),
    ('Pathogen',         str),
    ('Drug_class_CHN',   str),
))
# 毒力
sheet_VF = OrderedDict((
    ('Gene',             str),
    ('Reads_Num',        str),
    ('Coverage(%)',      str),
    ('Identity(%)',      str),
    ('Pathogen',         str),
    ('DatabaseID',       str),
    ('Annotation',       str),
))
# 物种鉴定
sheet_Pathogen = OrderedDict((
    ('Kingdom',         str),
    ('GenusName',       str),
    ('GenusCN',         str),
    ('GenusReads',      int),
    ('GroupName',       str),
    ('GroupReads',      int),
    ('ScientificName',  str),
    ('ChineseName',     str),
    ('Comment',         str),
    ('Reads_Number',    int),
    ('GMRN',            int),
    ('TPM',             float),
    ('RH-RPM',          float),
    ('Normalization',   int),
    ('Verification',    str),
    ('Focus',           str),
    ('Coverage',        float),
))
###################################################################################


def get_qc_df(in_dir, sample, library, rh_reads):
    # 获取质控前后数据的统计结果
    fastp_smry = Path(in_dir) / "01.QC/01.Fastp" / sample / library / f"{sample}.json"
    with open(fastp_smry, 'r') as fs_fp:
        fastp_stats = json.load(fs_fp)
    data = {
            'file': [f"{sample}_raw", f"{sample}_removeHost"],
            'Num reads': [fastp_stats['summary']['before_filtering']['total_reads'], rh_reads],
            'Num Bases': [fastp_stats['summary']['before_filtering']['total_bases'], rh_reads * 50]
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
        in_df = pd.read_excel(in_file, header=0)
        out_df = in_df[sheet_ARG.keys()]
    else:
        cols = list(sheet_ARG.keys())
        out_df = pd.DataFrame(columns=cols)
    return out_df


def get_vf_df(in_dir, sample, library):
    """ 毒力因子检测结果格式转换 """
    in_file = Path(in_dir) / "02.Annotation/03.VFDB" / sample / library / "VF.xls"
    if in_file.exists():
        out_df = pd.read_csv(in_file, sep='\t', header=0)
    else:
        cols = ['Gene', 'Reads_Num', 'Coverage(%)', 'Identity(%)', 'Pathogen', 'DatabaseID', 'Annotation']
        out_df = pd.DataFrame(columns=cols)
    return out_df


def get_tax_df(in_dir, sample):
    """ 病原鉴定结果格式转换 """
    tax_xlsx = Path(in_dir) / f'plugin_out/Result/{sample}/Total_Detail.xlsx'
    df = pd.read_excel(tax_xlsx)
    return df


@click.command(no_args_is_help=True)
@click.option('-w', '--workdir', help='Input dir, where the PIP outputs are.')
@click.option('-n', '--name', help='The sample name')
@click.option('-l', '--library', help='The library of the sample, DNA or RNA.')
def workflow(workdir, name, library):
    """ 将流程的分析结果按指定格式汇总为一个xls文件 """
    outdir = Path(workdir) / 'plugin_out/Result' / name
    if not Path(outdir).exists():
        os.system(f"mkdir -p {outdir}")

    # stat reads count of remove-host
    rh_file = Path(workdir) / '01.QC/02.RemoveHost' / name / library / f'{name}.fq.gz'
    if str(rh_file).endswith('gz'):
        with gzip.open(rh_file, 'rb') as rfp:
            n = 0
            for _ in rfp:
                n += 1
            rh_reads = int(n / 4)
    else:
        with open(rh_file, 'r') as rfp:
            n = 0
            for _ in rfp:
                n += 1
            rh_reads = int(n / 4)

    # merge sheet to xls
    sheet_name2df = OrderedDict()
    sheet_name2df['Data'], qc_stats = get_qc_df(workdir, name, library, rh_reads)
    sheet_name2df['ARG'] = get_arg_df(workdir, name, library)
    sheet_name2df['VF'] = get_vf_df(workdir, name, library)
    sheet_name2df['Pathogen'] = get_tax_df(workdir, name)

    out_xlsx = Path(outdir) / f"{name}.merge.xls"
    with pd.ExcelWriter(out_xlsx) as writer:
        for sheet_name, df in sheet_name2df.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    """ 根据fastp质控数据、去宿主后序列文件生成summary文件 """
    summary = defaultdict(list)
    summary['NumReads'].append(qc_stats['Num reads'][0])
    summary['NumBases'].append(qc_stats['Num Bases'][0])
    summary['NumSpcReads'].append(rh_reads)
    df_summary = pd.DataFrame(summary)
    out_summary = Path(outdir) / 'raw_summary.xls'
    df_summary.to_csv(out_summary, sep='\t', index=False)


if __name__ == '__main__':
    workflow()
