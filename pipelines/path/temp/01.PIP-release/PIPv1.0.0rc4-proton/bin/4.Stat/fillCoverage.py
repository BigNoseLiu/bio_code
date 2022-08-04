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

warnings.filterwarnings('ignore')

KINGDOM2PREFIX = {
    'bacteria': 'bacteria',
    'fungi': 'fungi',
    'LY': 'LY',
    'protozoa': 'protozoa',
    'viral': 'viral'
}


def fill_cov(xls, species2cov, key):
    """ 根据物种名获取其对应覆盖度 """
    df_xls = pd.read_csv(xls, sep='\t')
    if df_xls.empty:
        cols = list(df_xls.columns)
        cols.append('coverage')
        return pd.DataFrame(columns=cols)

    df_xls['coverage'] = df_xls.apply(lambda x: species2cov.get(x[key], 0), axis=1)
    if 'tax_id' in df_xls.columns:
        df_xls = df_xls.astype({'tax_id': int})
    return df_xls.fillna(0)


def copy_plotdir(source_dir, target_dir):
    """ 将目标目录拷贝至指定路径 """
    #if not Path(target_dir).exists():
    #    os.system(f'mkdir -p {str(target_dir)}')
    cmd = f'cp -r {str(source_dir)} {str(target_dir)}'
    if os.path.exists(source_dir):
        os.system(cmd)


def convert_str(data, tag):
    if data == '-':
        return '-'
    elif tag == 'float':
        return float(data)
    elif tag == 'int':
        return int(data)


def fill_genus_cname(df, db):
    # 补充检出物种对应的属中文名
    cols = df.columns

    name2genus_cn = {}
    with open(Path(db) / 'all.latin2chinese.INFO.txt', 'r') as rfp1:
        for _, line in enumerate(rfp1):
            ll = line.strip().split('\t')
            if ll:
                name2genus_cn[ll[2]] = ll[4]
    with open(Path(db) / 'GroupList.info', 'r') as rfp2:
        for _, line in enumerate(rfp2):
            ll = line.strip().split('\t')
            if ll:
                name2genus_cn[ll[1]] = ll[7]

    df['GenusCN'] = df.apply(lambda x: name2genus_cn.get(x['ScientificName']), axis=1)
    if 'GenusCN' not in cols:
        out_cols = cols.insert(2, 'GenusCN')
        return df[out_cols]
    else:
        return df[cols]


def fill_comment(df, src_file):
    cols = df.columns

    df_src = pd.read_csv(src_file, sep='\t')
    species2comment = {}
    for row in df_src[['ScientificName', 'Comment']].itertuples():
        species2comment[row[1]] = row[2]

    df['Comment'] = df.apply(lambda x: species2comment.get(x['ScientificName']), axis=1)

    if 'Comment' not in cols:
        out_cols = cols.insert(8, 'Comment')
        return df[out_cols]
    else:
        return df[cols]


@click.command(no_args_is_help=True)
@click.option('--db', help='Specified the DB version of Taxonomy')
@click.option('--indir', '-i', help='Input dir, where the PIP outputs are.')
@click.option('--name', '-n', help='The sample name')
@click.option('--library', '-l', help='The library of the sample, DNA or RNA.')
@click.option('--outdir', '-o', help='The output directory of results')
def main(db, indir, name, library, outdir):
    sample_dir = Path(indir) / '02.Annotation/01.Taxonomy' / name / library

    species2cov = defaultdict()
 
    for k, v in KINGDOM2PREFIX.items():
        cov = sample_dir / 'CovPlot' / f'{v}.{name}.cov'
        if cov.exists():
            df_cov = pd.read_csv(cov, sep='\t', header=None)
            if not df_cov.empty:
                df_cov.columns = ['Species', 'Component', 'Depth']
                df_sub = df_cov[df_cov['Depth'] > 0]
                df_gb = df_sub.groupby('Species').count()
                temp_sps2cov = defaultdict()
                for row in df_gb.itertuples():
                    temp_sps2cov[row[0]] = row[2]
            else:
                temp_sps2cov = defaultdict()
        else:
            temp_sps2cov = defaultdict()
        xls = sample_dir / f'{k}.xls'
        df_xls = fill_cov(xls, temp_sps2cov, 'species')
        df_xls.to_csv(xls, sep='\t', index=False)
        species2cov.update(temp_sps2cov)
        plot_dir = sample_dir / 'CovPlot' / f'{v}_readscount_coverage_plot'
        copy_plotdir(plot_dir, outdir)

    file_dict = {'Pathogeny_Detail.txt': 'Pathogeny_CommonSpecies.Reads.txt',
                 'Total_Detail.txt': 'Total_CommonSpecies.Reads.txt'}
    for file in file_dict.keys():
        file_obj = sample_dir / file
        # fill coverage
        df = fill_cov(file_obj, species2cov, 'ScientificName')
        if not df.empty:            
            # convert data type
            for col in ['GenusReads', 'GroupReads', 'Reads_Number', 'GMRN', 'Normalization']:
                df[col] = df.apply(lambda x: convert_str(x[col], 'int'), axis=1)
            for col in ['TPM', 'RH-RPM']:
                df[col] = df.apply(lambda x: convert_str(x[col], 'float'), axis=1)
            # fill genus chinese name
            df = fill_genus_cname(df, db)
            # fill comment
            src_file = Path(indir) / 'plugin_out/Result' / file_dict.get(file)
            df = fill_comment(df, src_file)
            df = df.rename(columns={'coverage': 'Coverage'})
        df.to_csv(file_obj, sep='\t', index=False)
        df.to_excel(file_obj.with_suffix('.xlsx'), index=False)


if __name__ == '__main__':
    main()
