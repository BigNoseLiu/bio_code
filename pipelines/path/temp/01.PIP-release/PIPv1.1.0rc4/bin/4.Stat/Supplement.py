#!/usr/bin/env python
import os
import sys
import json
import math
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
__DATE__ = "2022-06-14"
##########################################################################

KINGDOM2TAG = {
    'Bacteria': 'Bacteria',
    'Eukaryota:Fungi': 'Eukaryota_Fungi',
    'Eukaryota:Parasite': 'Eukaryota_Parasite',
    'Eukaryota:Protozoa': 'Eukaryota_Protozoa',
    'SpecialPathogens': 'SpecialPathogens',
    'Viruses': 'Viruses'
}


def create_dir(path):
    if not path.exists():
        os.system(f'mkdir -p {str(path)}')


def stat_common(df_sample, workdir, target_pathogen):
    """ 统计共有检出 """
    # 合并批次所有样本
    df_merge = pd.DataFrame()            # 未均一化结果
    df_merge_normal = pd.DataFrame()     # 已均一化结果
    spe2kingdom = {}
    spe2species_cn = {}
    for row in df_sample.itertuples():
        infile = Path(workdir) / 'plugin_out/Result/' / row[1] / 'Total_Detail.txt'
        if not infile.exists():
            continue
        df_in = pd.read_csv(infile, sep='\t', dtype=str)
        df_in.index = df_in['SpeciesSN']
        # 获取物种注释信息
        for line in df_in[['Kingdom', 'SpeciesCN']].itertuples():
            spe2kingdom[line[0]] = line[1]
            spe2species_cn[line[0]] = line[2]
        # 获取物种检出reads及合并数据
        df_sub = df_in[['SpeciesReads']].rename(columns={'SpeciesReads': row[1]}).astype(int)
        df_sub_normal = df_in[['SpeciesReads(Normalization)']] \
            .rename(columns={'SpeciesReads(Normalization)': row[1]}).astype(int)
        df_merge = pd.concat([df_merge, df_sub], axis=1)
        df_merge_normal = pd.concat([df_merge_normal, df_sub_normal], axis=1)
    df_merge = df_merge.fillna(0)
    df_merge_normal = df_merge_normal.fillna(0)

    # 统计批次中物种检出频率及排位
    idx2smp = {i: col for i, col in enumerate(df_merge.columns)}
    spe2smp2comment = defaultdict(defaultdict)
    spe2comment = {}
    total_smp = len(df_sample)
    for row in df_merge.itertuples():
        # 统计检出频率
        common = f'{str(len([i for i in row[1:] if i != 0]))}/{total_smp}'
        # 统计排位
        value2rank = {}
        for i, v in enumerate(sorted(set(row[1:]), reverse=True), 1):
            value2rank[v] = i
        for i, v in enumerate(row[1:]):
            spe2smp2comment[row[0]][idx2smp.get(i)] = f'{common}|{str(value2rank.get(v))}'
            spe2comment[row[0]] = f'{common}|{str(value2rank.get(v))}'

    # 补充物种注释信息
    out_cols = ['Kingdom', 'SpeciesSN', 'SpeciesCN', 'Comment'] + list(df_merge.columns)
    for df in df_merge, df_merge_normal:
        df['Kingdom'] = [spe2kingdom.get(i, '-') for i in df.index]
        df['SpeciesSN'] = df.index
        df['SpeciesCN'] = [spe2species_cn.get(i, '-') for i in df.index]
        df['Comment'] = [spe2comment.get(i, '-').split('|')[0] for i in df.index]

    outfile = Path(workdir) / 'plugin_out/Result' / 'Total_CommonSpecies.Reads.txt'
    df_merge.to_csv(outfile, sep='\t', index=False, columns=out_cols)
    with pd.ExcelWriter(outfile.with_suffix('.xlsx')) as writer1:
        df_merge.to_excel(writer1, sheet_name='批次物种检出表', index=False, columns=out_cols)
        df_ptg = df_merge[df_merge['SpeciesSN'].isin(target_pathogen)]
        df_ptg.to_excel(writer1, sheet_name='批次常见菌检出矩阵表', index=False, columns=out_cols)

    outfile_nor = Path(workdir) / 'plugin_out/Result' / 'Total_CommonSpecies.Normalize.txt'
    df_merge_normal.to_csv(outfile_nor, sep='\t', index=False, columns=out_cols)
    with pd.ExcelWriter(outfile_nor.with_suffix('.xlsx')) as writer2:
        df_merge_normal.to_excel(writer2, sheet_name='批次物种检出表', index=False, columns=out_cols)
        df_ptg_nor = df_merge_normal[df_merge_normal['SpeciesSN'].isin(target_pathogen)]
        df_ptg_nor.to_excel(writer2, sheet_name='批次常见菌检出矩阵表', index=False, columns=out_cols)

    return spe2comment


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


def get_nc(df_samples, n_million, workdir):
    """ 获取阴性质控样本检出结果 """
    nc_smp = ''
    for row in df_samples.itertuples():
        if row[1].startswith('Dzk') and 'DH' not in row[1]:
            nc_smp = row[1]
            nc_lib = row[2]
    if not nc_smp:
        nc_smp = df_samples.iloc[0, 0]
        nc_lib = df_samples.iloc[0, 1]
        print(f"Warning: NC sample name's format is unvalidated, replace with {nc_smp}")

    qc_file = Path(workdir) / '01.QC/01.Fastp' / nc_smp / nc_lib / f'{nc_smp}.json'
    with open(qc_file, 'r') as rfp:
        qc_info = json.load(rfp)
    coe = int(n_million) * 10 ** 6 / qc_info["summary"]["before_filtering"]["total_reads"]

    nc_file = Path(workdir) / 'plugin_out/Result' / nc_smp / 'Total_Detail.txt'
    df_nc = pd.read_csv(nc_file, sep='\t', dtype=str)
    df_nc['SpeciesReads(Normalization)'] = [math.ceil(int(i) * coe) for i in df_nc['SpeciesReads']]
    nc_spe2rc = {row[1]: int(row[2]) for row in df_nc[['SpeciesSN', 'SpeciesReads(Normalization)']].itertuples()}

    return nc_spe2rc


def fill_normalization(df_data, n_million, qc_file):
    """ 计算均一化reads并输出至结果表 """
    with open(qc_file, 'r') as rfp:
        qc_info = json.load(rfp)
    coe = int(n_million) * 10 ** 6 / qc_info["summary"]["before_filtering"]["total_reads"]

    df_data['GenusReads(Normalization)'] = [math.ceil(int(i) * coe) for i in df_data['GenusReads']]
    df_data['GroupReads(Normalization)'] = [math.ceil(int(i) * coe) for i in df_data['GroupReads']]
    df_data['SpeciesReads(Normalization)'] = [math.ceil(int(i) * coe) for i in df_data['SpeciesReads']]

    return df_data


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


def fill_common(file_path, spe2comment):
    """ 补充共有检出统计 """
    df_data = pd.read_csv(file_path, sep='\t', dtype=str)

    if df_data.empty:
        print(f'Warning: {infile} is empty, please check it.')
        return

    if not spe2comment:
        df_data['Comment'] = ['-' for i in range(len(df_data))]
    else:
        df_data['Comment'] = [spe2comment.get(i, '0') for i in df_data['SpeciesSN']]

    file_lst = [str(file_path), str(file_path.with_suffix('.xlsx'))]
    os.system('rm ' + ' '.join(file_lst))
    df_data.to_csv(file_path, sep='\t', index=False)
    df_data.to_excel(file_path.with_suffix('.xlsx'), index=False)

    return


def fill(file_path, normalize, qc_file, cov_file, bg_info, nc_spe2rc):
    """ 补充背景菌、覆盖度、阴控数据 """
    if not file_path.parent.exists():
        os.system(f'mkdir -p {str(outfile.parent)}')

    # Data Input
    df_data = pd.read_csv(file_path, sep='\t', dtype=str)
    if df_data.empty:
        print(f'Warning: {infile} is empty, please check it.')
        return

    # Fill normalization
    if qc_file.exists():
        df_data = fill_normalization(df_data, normalize, qc_file)
    else:
        print(f'Warning: {qc_file} is not exists, please check it.')

    # Fill background species stat
    if len(bg_info):
        df_data = fill_background(df_data, bg_info)
    else:
        print(f'Warning: Background Stat Data is empty, please check it.')

    # Fill coverage
    if cov_file.exists():
        with open(cov_file, 'r') as rfp:
            spe2cov = json.load(rfp)
        df_data = fill_cov(df_data, spe2cov)
    else:
        print(f'Warning: {cov_file} not exists, please check it.')

    # Fill negative control sample's data
    if len(nc_spe2rc):
        df_data = fill_nc(df_data, nc_spe2rc)
    else:
        print('Warning: Negative Control Data is empty, please check it.')

    tmp_csv = file_path.with_suffix('.tmp')
    df_data.to_csv(tmp_csv, sep='\t', index=False)
    os.system(f'mv {str(tmp_csv)} {str(file_path)}')
    df_data.to_excel(file_path.with_suffix('.xlsx'), index=False)


@click.command(no_args_is_help=True)
@click.option('-w', '--workdir', help='Where the PIP output dir is <str>')
@click.option('-s', '--samples', help='Sample list file <str>')
@click.option('--nc', help='Negative Control sample\'s name <str>', default='')
@click.option('--db', help='PIDB path <str>')
@click.option('-b', '--bg_ver', help='Background database\'s version. <str>', default='BG.v1')
@click.option('-n', '--normalize', help='n M reads to normalize. [20]', default=20)
@click.option('-t', '--thread', help='How many threads can be used. [21]', default=21)
def main(workdir, samples, nc, db, bg_ver, normalize, thread):
    """ 补充阴控、背景菌等相关数据 """
    # 样本列表
    df_smp = pd.read_csv(samples, sep='\t', dtype=str, header=None)

    # 获取阴性样本数据
    nc_spe2rc = get_nc(df_smp, normalize, workdir)

    # 获取背景菌数据库数据
    bg_info = {}
    pool_bg = Pool(int(thread))
    for smp_site in set(df_smp[len(df_smp.columns) - 1]):
        bg_file = Path(db) / bg_ver / '.'.join([smp_site.lower(), 'xlsx'])
        bg_info[smp_site] = pool_bg.apply_async(get_background_info, (bg_file, )).get()
    pool_bg.close()
    pool_bg.join()

    qc_dir = Path(workdir) / '01.QC/01.Fastp'
    plg_dir = Path(workdir) / 'plugin_out/Result'

    # 补充覆盖度、背景菌数据、共有检出等至结果表
    pool_fill = Pool(int(thread))
    for file in 'Total_Detail.txt', 'Pathogeny_Detail.txt':
        for row in df_smp.itertuples():
            tax_file = plg_dir / row[1] / file
            if not tax_file.exists():
                print(f'Error: {tax_file} not exists, please check it.')
                continue
            qc_file = qc_dir / row[1] / row[2] / f'{row[1]}.json'
            sub_bg_info = bg_info[row[-1]]
            cov_file = Path(workdir) / f'02.Annotation/01.Taxonomy/{row[1]}/{row[2]}/Coverage/03.stat/stat.json'
            pool_fill.apply_async(fill, args=(tax_file, normalize, qc_file, cov_file, sub_bg_info, nc_spe2rc))
    pool_fill.close()
    pool_fill.join()

    # 获取3, 4级致病物种列表
    pc_file = Path(db) / 'PathogenicClassificationList.txt'
    df_pc = pd.read_csv(pc_file, sep='\t', dtype=str, header=1)
    df_pc = df_pc[(df_pc['致病综合打分'] == '3') | (df_pc['致病综合打分'] == '4')]
    target_pathogen = set(df_pc['l_name_s'])

    # 批次统计：共有物种检出
    spe2comment = stat_common(df_smp, workdir, target_pathogen)
    if len(spe2comment) == 0:
        print('Warning: Common Stat Data is empty, please check it.')

    # Fill comment of common species
    for file in 'Total_Detail.txt', 'Pathogeny_Detail.txt':
        for row in df_smp.itertuples():
            tax_file = plg_dir / row[1] / file
            if not tax_file.exists():
                print(f'Error: {str(tax_file)} not exists, please check it.')
                continue
            else:
                fill_common(tax_file, spe2comment)


if __name__ == '__main__':
    main()
