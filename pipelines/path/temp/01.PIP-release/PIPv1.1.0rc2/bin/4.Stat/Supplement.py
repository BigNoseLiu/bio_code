#!/usr/bin/env python
import os
import sys
import json
import click
import warnings
import subprocess
import pandas as pd
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../conf'))
import configure as cfg

warnings.filterwarnings('ignore')

############################ developer info ##############################
__AUTHOR__ = "Sujiawei"
__EMAIL__ = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__ = "2022-05-31"
##########################################################################

KINGDOM2TAG = {
    'Bacteria': 'bacteria',
    'Eukaryota:Fungi': 'fungi',
    'Eukaryota:Parasite': 'parasite',
    'Eukaryota:Protozoa': 'protozoa',
    'SpecialPathogens': 'special_pathogen',
    'Viruses': 'virus'
}


def create_dir(path):
    if not path.exists():
        os.system(f'mkdir -p {str(path)}')


def get_background_info(bg_file):
    """ 获取背景数据库相关数据 """
    df_bg_n = pd.read_excel(bg_file, dtype=str, sheet_name='csf20M均一化')
    df_bg_g = pd.read_excel(bg_file, dtype=str, sheet_name='csf属排名')
    df_bg_s = pd.read_excel(bg_file, dtype=str, sheet_name='csf种排名')
    gns2bg_rank = {}
    spe2bg_rank = {}
    spe2bg_lvl = {}
    spe2bg_fre = {}
    for row in df_bg_n[['ScientificName', 'Freq/[Min-25%-Median-75%-Max](Historically)',
                        'Background_level']].itertuples():
        spe2bg_fre[row[1]] = row[2]
        spe2bg_lvl[row[1]] = row[3]
    for row in df_bg_g[['Genus', '[Min-25%-Median-75%-Max](Historically)']].itertuples():
        gns2bg_rank[row[1]] = row[2]
    for row in df_bg_s[['ScientificName', '[Min-25%-Median-75%-Max](Historically)']].itertuples():
        spe2bg_rank[row[1]] = row[2]

    return {'gns2bg_rank': gns2bg_rank,
            'spe2bg_rank': spe2bg_rank,
            'spe2bg_lvl': spe2bg_lvl,
            'spe2bg_fre': spe2bg_fre}


def get_nc(df_samples, workdir):
    """ 获取阴性质控样本检出结果 """
    nc_smp = ''
    for row in df_samples.itertuples():
        if row[1].startswith('Dzk') and 'DH' not in row[1]:
            nc_smp = row[1]
    if not nc_smp:
        nc_smp = df_samples.iloc[0, 0]
        print(f"NC sample name's format is unvalidated, replace with {nc_smp}")

    nc_file = Path(workdir) / 'plugin_out/Result' / nc_smp / 'Total_Detail.txt'
    df_nc = pd.read_csv(nc_file, sep='\t', dtype=str)
    nc_spe2rc = {row[1]: int(row[2]) for row in df_nc[['SpeciesSN', 'SpeciesReads(Normalization)']].itertuples()}
    return nc_spe2rc


def fill_background(df_data, bg_info):
    """ 补充背景菌数据 """
    df_data['Freq/[Min-25%-Median-75%-Max](Historically)'] = [bg_info['spe2bg_fre'].get(i, '-') for i in df_data['SpeciesSN']]
    df_data['Background_level'] = [bg_info['spe2bg_lvl'].get(i, '-') for i in df_data['SpeciesSN']]

    genus_rank = []
    for row in df_data[['GenusSN', 'GenusRank']].itertuples():
        new_rank = '|'.join([row[2], bg_info['gns2bg_rank'].get(row[1], '-')])
        genus_rank.append(new_rank)
    df_data['GenusRank'] = genus_rank

    spe2rank = {}
    for kd in KINGDOM2TAG.keys():
        df_temp = df_data[df_data['Kingdom'] == kd].astype({'SpeciesReads': 'int'})
        values = sorted(set(df_temp['SpeciesReads']), reverse=True)
        value2rank = {v: i for i, v in enumerate(values, 1)}
        for row in df_temp[['SpeciesSN', 'SpeciesReads']].itertuples():
            spe2rank[row[1]] = '|'.join([str(value2rank.get(row[2])), bg_info['spe2bg_rank'].get(row[1], '-')])

    df_data['SpeciesRank'] = [spe2rank.get(i) for i in df_data['SpeciesSN']]

    return df_data


def fill_cov(df, spe2cov):
    """ 补充覆盖度 """
    for col in 'Coverage', 'BinCoverage':
        value_lst = []
        for species in df['SpeciesSN']:
            if spe2cov.get(species):
                value_lst.append(spe2cov[species][col])
            else:
                value_lst.append(0)
        df[col] = value_lst
    return df


def fill_nc(df_data, nc_spe2rc):
    """ 补充阴性质控数据 """
    if not nc_spe2rc:
        return df_data

    df_data['NCCount'] = [nc_spe2rc.get(i, 0) for i in df_data['SpeciesSN']]

    nc_fdc = []
    for row in df_data[['NCCount', 'SpeciesReads(Normalization)']].itertuples():
        if row[1] == 0:
            nc_fdc.append('0.000')
        else:
            nc_fdc.append('{:.3f}'.format(int(row[2]) / row[1]))
    df_data['FoldChange(NC)'] = nc_fdc

    return df_data


def fill(file_path, bg_info, cov_file, nc_spe2rc):
    """ 补充背景菌、覆盖度、阴控数据 """
    if not file_path.parent.exists():
        os.system(f'mkdir -p {str(outfile.parent)}')

    df_file = pd.read_csv(file_path, sep='\t', dtype=str)
    if df_file.empty:
        print(f'{infile} is empty, please check it.')
        return
    df_data = fill_background(df_file, bg_info)

    if cov_file.exists():
        with open(cov_file, 'r') as rfp:
            spe2cov = json.load(rfp)
        df_data = fill_cov(df_data, spe2cov)
    else:
        print(f'{cov_file} not exists, please check it.')

    if len(nc_spe2rc):
        df_data = fill_nc(df_data, nc_spe2rc)
    else:
        print('Negative Control Data is empty, please check it.')

    tmp_csv = file_path.with_suffix('.tmp')
    df_data.to_csv(tmp_csv, sep='\t', index=False)
    os.system(f'mv {str(tmp_csv)} {str(file_path)}')
    df_data.to_excel(file_path.with_suffix('.xlsx'), index=False)


def stat_common(df_sample, workdir):
    """ 统计共有检出 """
    df_merge = pd.DataFrame()
    for row in df_sample.itertuples():
        infile = Path(workdir) / 'plugin_out/Result/' / row[1] / 'Total_Detail.txt'
        df_in = pd.read_csv(infile, sep='\t', dtype=str)
        df_in.index = df_in['SpeciesSN']
        df_sub = df_in[['SpeciesReads']].rename(columns={'SpeciesReads': row[1]}).astype(int)
        df_merge = pd.concat([df_merge, df_sub], axis=1)
    df_merge = df_merge.fillna(0)

    data = []
    for row in df_merge.itertuples():
        common = f'{str(len([i for i in row[1:] if i != 0]))}/{str(len(row) - 1)}'
        data_row = [row[0], common]
        value2rank = {}
        for i, v in enumerate(sorted(set(row[1:]), reverse=True), 1):
            value2rank[v] = i
        for v in row[1:]:
            data_row.append(value2rank.get(v))
        data.append(data_row)
    data_cols = ['SpeciesSN', 'Common'] + list(df_sample[0])
    df_common = pd.DataFrame(data, columns=data_cols)

    outfile = Path(workdir) / 'plugin_out/Result' / 'Total_CommonSpecies.Stat.txt'
    df_common.to_csv(outfile, sep='\t', index=False)
    df_common.to_excel(outfile.with_suffix('.xlsx'), index=False)


@click.command(no_args_is_help=True)
@click.option('-w', '--workdir', help='Where the PIP output dir is <str>')
@click.option('-s', '--samples', help='Sample list file <str>')
@click.option('--nc', help='Negative Control sample\'s name <str>', default='')
@click.option('--db', help='PIDB path <str>')
@click.option('-b', '--bg_ver', help='Background database\'s version. <str>', default='BG.v1')
@click.option('-t', '--thread', help='How many threads can be used. [21]', default=21)
def main(workdir, samples, nc, db, bg_ver, thread):
    """ 补充阴控、背景菌等相关数据 """
    # 样本列表
    df_smp = pd.read_csv(samples, sep='\t', dtype=str, header=None)

    # 获取阴性样本数据
    nc_spe2rc = get_nc(df_smp, workdir)

    # 获取背景菌数据库数据
    bg_info = {}
    pool_bg = Pool(int(thread))
    for smp_site in set(df_smp[len(df_smp.columns) - 1]):
        bg_file = Path(db) / bg_ver / '.'.join([smp_site.lower(), 'xlsx'])
        bg_info[smp_site] = pool_bg.apply_async(get_background_info, (bg_file, )).get()
    pool_bg.close()
    pool_bg.join()

    plg_dir = Path(workdir) / 'plugin_out/Result'

    # 补充覆盖度、背景菌数据、共有检出等至结果表
    pool_fill = Pool(int(thread))
    for file in 'Total_Detail.txt', 'Pathogeny_Detail.txt':
        for row in df_smp.itertuples():
            file_path = plg_dir / row[1] / file
            cov_file = Path(workdir) / f'02.Annotation/01.Taxonomy/{row[1]}/{row[2]}/Coverage/03.stat/stat.json'
            sub_bg_info = bg_info[row[-1]]
            pool_fill.apply_async(fill, args=(file_path, sub_bg_info, cov_file, nc_spe2rc))
    pool_fill.close()
    pool_fill.join()

    # 批次统计：共有物种检出
    stat_common(df_smp, workdir)


if __name__ == '__main__':
    main()