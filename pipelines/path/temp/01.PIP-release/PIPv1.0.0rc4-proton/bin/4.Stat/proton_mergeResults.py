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
# 物种鉴定-非病毒
tax_out_cols_dtypes = OrderedDict((
    ('genus reads cuont',      str),
    ('coverage',               str),
    ('maxscore',               str),
    ('specise',                str),
    ('reads count',            str)
))
# 物种鉴定-病毒
viral_out_cols_dtypes = OrderedDict((
    ('coverage',               str),
    ('maxscore',               str),
    ('specise',                str),
    ('reads count',            str)
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


def get_tax_df(in_dir, sample, library):
    """ 病原鉴定结果格式转换 """
    # 定义输出结果的列名
    ctg2cols = OrderedDict((
        ('bacteria.xls',    list(tax_out_cols_dtypes.keys()) + ['', '', '', '', '', '']),
        ('fungi.xls',       list(tax_out_cols_dtypes.keys()) + ['', '', '', '', '', '']),
        ('viral.xls',       list(viral_out_cols_dtypes.keys()) + ['', '', '', '', '', '']),
        ('LY.xls',          list(tax_out_cols_dtypes.keys()) + ['', '', '', '', '', '']),
        ('protozoa.xls',    list(tax_out_cols_dtypes.keys()) + ['', '', '', '', '', ''])
    ))
    # 定义输出结果的表名
    xls2sheetname = OrderedDict((
        ('bacteria.xls',    'bacteria2.pick'),
        ('fungi.xls',       'fungi2.pick'),
        ('viral.xls',       'viral2.pick50'),
        ('LY.xls',          'LY2.pick'),
        ('protozoa.xls',    'protozoa2.pick')
    ))

    category2df = OrderedDict()
    # 输入目录
    xls_dir = Path(in_dir) / "02.Annotation/01.Taxonomy" / sample / library
    for xls, cols in ctg2cols.items():
        infile = xls_dir / xls
        df = pd.read_csv(infile, sep='\t', header=0, dtype=str)
        # 获取输入文件的指定列
        sub_df = df.iloc[:, [3, 9, 1, 2, 5, 4, 6, 8, 7, 10]]
        if xls != 'viral.xls':
            sub_df.insert(0, 'genus reads cuont', list('0' * len(sub_df)))
            sub_df['genus reads cuont'] = sub_df['genus reads cuont'].astype(int)
        try:
            sub_df.columns = cols
        except:
            print(sub_df.columns)
            print(xls)
            sys.exit(1)
        sub_df['reads count'] = sub_df['reads count'].astype(int)

        # 按输出格式修改列名
        sheet_name = xls2sheetname[xls]
        category2df[sheet_name] = sub_df

    return category2df


def generate_summary(qc_stats, tax_dfs):
    """ 根据fastp质控数据、物种鉴定结果生成summary文件 """
    summary = defaultdict(list)
    summary['NumReads'].append(qc_stats['Num reads'][0])
    summary['NumBases'].append(qc_stats['Num Bases'][0])
    species_reads_cnt = 0
    for category, df in tax_dfs.items():
        cnt_temp = sum(df['reads count'].astype(int))
        species_reads_cnt += cnt_temp
    summary['NumSpcReads'].append(species_reads_cnt)

    df_summary = pd.DataFrame(summary)
    return df_summary


@click.command(no_args_is_help=True)
@click.option('--indir', '-i', help='Input dir, where the PIP outputs are.')
@click.option('--name', '-n', help='The sample name')
@click.option('--library', '-l', help='The library of the sample, DNA or RNA.')
@click.option('--outdir', '-o', help='The output directory of merged result of summaries')
def workflow(indir, name, library, outdir):
    """ 将流程的分析结果按指定格式汇总为一个xls文件 """
    if not Path(outdir).exists():
        os.system(f"mkdir -p {outdir}")

    # merge sheet to xls
    sheet_name2df = OrderedDict()
    #sheet_name2df['data'], qc_stats = get_qc_df(indir, name, library) #hy
    sheet_name2df['ARG'] = get_arg_df(indir, name, library)
    sheet_name2df['VF'] = get_vf_df(indir, name, library)
    tax_dfs = get_tax_df(indir, name, library)
    sheet_name2df.update(tax_dfs)
    with pd.ExcelWriter(Path(outdir) / f"{name}.merge.xls") as writer:
        for sheet_name, df in sheet_name2df.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    #df_summary = generate_summary(qc_stats, tax_dfs) #hy
    #out_summary = Path(outdir) / 'raw_summary.xls' #hy
    #df_summary.to_csv(out_summary, sep='\t', index=False) #hy


if __name__ == '__main__':
    workflow()
