#!/usr/bin/env python
import os
import re
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
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as pat

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


def get_quotient(numerator, denominator):
    if denominator == 0:
        return 0
    else:
        return numerator / denominator


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
    spe2smp2ratio = defaultdict(dict)
    spe2smp2comment = defaultdict(dict)
    spe2comment = {}
    total_smp = len(df_sample)

    for row in df_merge.itertuples():
        # 统计样本物种检出与批次物种检出最大值的比值
        batch_max = max(row[1:])
        for i, v in enumerate(row[1:]):
            spe2smp2ratio[row[0]][idx2smp.get(i)] = '{:.2f}'.format(get_quotient(v, batch_max))

    for row in df_merge.itertuples():
        # 统计检出频率
        common = f'{str(len([i for i in row[1:] if i != 0]))}/{total_smp}'
        spe2comment[row[0]] = common
        # 统计排位
        value2rank = {}
        for i, v in enumerate(sorted(set(row[1:]), reverse=True), 1):
            value2rank[v] = i
        for i, v in enumerate(row[1:]):
            smp = idx2smp.get(i)
            spe2smp2comment[row[0]][smp] = f'{common}|{str(value2rank.get(v))}'
            ratio = spe2smp2ratio[row[0]].get(smp)
            spe2smp2comment[row[0]][smp] = f'{common}|{str(value2rank.get(v))}({ratio})'

    # 补充物种注释信息
    out_cols = ['Kingdom', 'SpeciesSN', 'SpeciesCN', 'Comment'] + list(df_merge.columns)
    for df in df_merge, df_merge_normal:
        df['Kingdom'] = [spe2kingdom.get(i, '-') for i in df.index]
        df['SpeciesSN'] = df.index
        df['SpeciesCN'] = [spe2species_cn.get(i, '-') for i in df.index]
        df['Comment'] = [spe2comment.get(i, '-') for i in df.index]

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

    return spe2smp2comment, df_merge, spe2kingdom


def get_background_info(bg_file):
    """ 获取背景数据库相关数据 """
    df_bg_n = pd.read_excel(bg_file, dtype=str, sheet_name='20M均一化')
    df_bg_g = pd.read_excel(bg_file, dtype=str, sheet_name='属排名')
    df_bg_s = pd.read_excel(bg_file, dtype=str, sheet_name='种排名')
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
    tag2smp_lib = {}
    for row in df_samples.itertuples():
        if row[1].startswith('Dzk') and 'DH' not in row[1]:
            # DNA 阴性质控样本
            tag2smp_lib['Dzk'] = '_'.join([row[1], row[2]])
        elif row[1].startswith('Rzk') and 'DH' not in row[1]:
            # RNA 阴性质控样本
            tag2smp_lib['Rzk'] = '_'.join([row[1], row[2]])
        elif row[1].startswith('Dzk') and 'DH' in row[1]:
            # DNA 阴性质控样本-去宿主处理
            tag2smp_lib['Dzk_DH'] = '_'.join([row[1], row[2]])
        elif row[1].startswith('Rzk') and 'DH' in row[1]:
            # RNA 阴性质控样本-去宿主处理
            tag2smp_lib['Rzk_DH'] = '_'.join([row[1], row[2]])

    if not tag2smp_lib:
        smp_lib = '_'.join([df_samples.iloc[0, 0], df_samples.iloc[0, 1]])
        tag2smp_lib['Dzk'] = smp_lib
        tag2smp_lib['Rzk'] = smp_lib
        tag2smp_lib['Dzk_DH'] = smp_lib
        tag2smp_lib['Rzk_DH'] = smp_lib
        print(f"Warning: NC sample name's format is unvalidated, replace with {smp_lib}")

    nc_spe2tag2rc = defaultdict(dict)
    for tag, smp_lib in tag2smp_lib.items():
        smp = smp_lib.split('_')[0]
        lib = smp_lib.split('_')[1]
        # get coefficient
        qc_file = Path(workdir) / '01.QC/01.Fastp' / smp / lib / f'{smp}.json'
        with open(qc_file, 'r') as rfp:
            qc_info = json.load(rfp)
        coe = get_quotient(int(n_million) * 10 ** 6, qc_info["summary"]["before_filtering"]["total_reads"])
        nc_file = Path(workdir) / 'plugin_out/Result' / smp / 'Total_Detail.txt'
        df_nc = pd.read_csv(nc_file, sep='\t', dtype=str)
        df_nc['SpeciesReads(Normalization)'] = [math.ceil(int(i) * coe) for i in df_nc['SpeciesReads']]
        # get normalized reads
        for row in df_nc[['SpeciesSN', 'SpeciesReads(Normalization)']].itertuples():
            nc_spe2tag2rc[row[1]][tag] = int(row[2])

    return nc_spe2tag2rc


def sort_filter(df_data, target_pathogen, top=10):
    """ 按界进行归类，再按属reads、种reads分别由高到低进行排序; 根据致病菌分级列表进行筛选, 每个属最多保留前top个物种结果 """
    df_ret = pd.DataFrame(columns=df_data.columns)
    for kd in KINGDOM2TAG.keys():
        df_kd = df_data[df_data['Kingdom'] == kd]
        if df_kd.empty:
            continue
        df_kd = df_kd.astype({'GenusReads': int, 'SpeciesReads': int})
        df_sort = df_kd.sort_values(by=['GenusReads', 'SpeciesReads'], ascending=False)
        df_filter = pd.DataFrame(columns=df_data.columns)
        for gns in set(df_sort['GenusSN']):
            df_gns = df_sort[df_sort['GenusSN'] == gns]
            # is pathogens
            df_is_ptg = df_gns[df_gns['SpeciesSN'].isin(target_pathogen)]
            df_filter = pd.concat([df_filter, df_is_ptg])
            if len(df_is_ptg) >= 10:
                continue
            # is not pathogens
            df_not_ptg = df_gns[~df_gns['SpeciesSN'].isin(target_pathogen)]
            df_filter = pd.concat([df_filter, df_not_ptg.iloc[:top-len(df_is_ptg)]])
        df_filter = df_filter.sort_values(by=['GenusReads', 'SpeciesReads'], ascending=False)
        df_ret = pd.concat([df_ret, df_filter])
    return df_ret


def fill_db_code(df_data, spe2db_code):
    """ 补充DatabaseCode """
    df_data['DatabaseCode'] = [spe2db_code.get(i, '-') for i in df_data['SpeciesSN']]
    return df_data


def fill_normalization(df_data, n_million, qc_file):
    """ 计算均一化reads并输出至结果表 """
    with open(qc_file, 'r') as rfp:
        qc_info = json.load(rfp)
    coe = get_quotient(int(n_million) * 10 ** 6, qc_info["summary"]["before_filtering"]["total_reads"])

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


def fill_nc(smp, lib, df_data, nc_spe2tag2rc):
    """ 补充阴性质控数据 """
    if not nc_spe2tag2rc:
        return df_data

    if 'DH' not in smp:
        if lib == 'DNA':
            tag = 'Dzk'
        elif lib == 'RNA':
            tag = 'Rzk'
    elif 'DH' in smp:
        if lib == 'DNA':
            tag = 'Dzk_DH'
        elif lib == 'RNA':
            tag = 'Rzk_DH'
    else:
        tag = 'Dzk'

    nc_counts = []
    for spe in df_data['SpeciesSN']:
        if spe in nc_spe2tag2rc:
            nc_counts.append(nc_spe2tag2rc[spe].get(tag, 0))
        else:
            nc_counts.append(0)
    df_data['NCCount'] = nc_counts

    nc_fdc = []
    for row in df_data[['NCCount', 'SpeciesReads(Normalization)']].itertuples():
        if row[1] == 0:
            nc_fdc.append('0.000')
        else:
            nc_fdc.append('{:.3f}'.format(int(get_quotient(row[2], row[1]))))
    df_data['FoldChange(NC)'] = nc_fdc

    return df_data


def fill_common(file_path, sample, spe2smp2comment):
    """ 补充共有检出统计 """
    df_data = pd.read_csv(file_path, sep='\t', dtype=str)

    if df_data.empty:
        print(f'Warning: {infile} is empty, please check it.')
        return

    if not spe2smp2comment:
        df_data['Comment'] = ['-' for i in range(len(df_data))]
    else:
        comment = []
        for spe in df_data['SpeciesSN']:
            if spe in spe2smp2comment:
                comment.append(spe2smp2comment[spe].get(sample, '-'))
            else:
                comment.append('-')
        df_data['Comment'] = comment

    file_lst = [str(file_path), str(file_path.with_suffix('.xlsx'))]
    os.system('rm ' + ' '.join(file_lst))
    df_data.to_csv(file_path, sep='\t', index=False)
    df_data.to_excel(file_path.with_suffix('.xlsx'), index=False)

    return


def fill(file_path, normalize, qc_file, cov_file, bg_info, nc_spe2tag2rc, spe2db_code, smp, lib, target_pathogen):
    """ 补充背景菌、覆盖度、阴控数据 """
    if not file_path.parent.exists():
        os.system(f'mkdir -p {str(outfile.parent)}')

    # Data Input
    df_data = pd.read_csv(file_path, sep='\t', dtype=str)
    if df_data.empty:
        print(f'Warning: {infile} is empty, please check it.')
        return

    # Sort and filter data
    df_data = sort_filter(df_data, target_pathogen)

    # Fill DatabaseCode
    if len(spe2db_code):
        df_data = fill_db_code(df_data, spe2db_code)
    else:
        print(f'Warning: DatabaseCode is empty, please check it.')

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
    if len(nc_spe2tag2rc):
        df_data = fill_nc(smp, lib, df_data, nc_spe2tag2rc)
    else:
        print('Warning: Negative Control Data is empty, please check it.')

    tmp_csv = file_path.with_suffix('.tmp')
    df_data.to_csv(tmp_csv, sep='\t', index=False)
    os.system(f'mv {str(tmp_csv)} {str(file_path)}')
    df_data.to_excel(file_path.with_suffix('.xlsx'), index=False)


def draw_plot(x_data, y_data, name, sample, outfile):
    """ 绘制批次强阳物种统计柱状图 """
    fig, ax = plt.subplots(figsize=(12, 12))
    color_lst = []
    smp_idx = 0
    for i, smp in enumerate(x_data):
        if smp == sample:
            smp_idx = i
            color_lst.append('coral')
        else:
            color_lst.append('dodgerblue')
    ax.barh(x_data, y_data, color=color_lst, )

    # 坐标轴设置
    ax.set_xlabel('Unique Mapping Reads', fontdict={'fontsize': 14})
    ax.set_ylabel('Sample ID', fontdict={'fontsize': 14})
    ax.set_xlim(0, max(y_data))

    # 处理边距
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('dodgerblue')
    ax.spines['bottom'].set_position(('data', 0))
    ax.margins(False)

    # 处理图例、标题
    plt.subplots_adjust(top=0.77, bottom=0.14)
    plt.suptitle(f"{name}", fontsize=19, fontweight='semibold', y=0.95)
    patches = [pat.Patch(color='dodgerblue', label='Unique Mapping Reads'),
               pat.Patch(color='coral', label='Current Sample ID')]
    plt.gca().legend(handles=patches, loc='upper center', ncol=2, fontsize=12, bbox_to_anchor=(0.5, 1.19))

    # 调整图片位置
    plt.gcf().subplots_adjust(right=0.875, bottom=0.2)

    # 坐标轴刻度设置
    if smp_idx:
        ax.get_yticklabels()[smp_idx].set_color('coral')

    plt.savefig(outfile)


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
    nc_spe2tag2rc = get_nc(df_smp, normalize, workdir)

    # 获取Database code
    spe2db_code = {}
    with open(Path(db) / 'all.Taxonomy.txt') as rfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                if ll[4] == 'Viruses':
                    spe2db_code[ll[12]] = ll[0]
                else:
                    spe2db_code[ll[10]] = ll[0]

    # 获取背景菌数据库数据
    bg_info = {}
    pool_bg = Pool(int(thread))
    for smp_site in set(df_smp[len(df_smp.columns) - 1]):
        bg_file = Path(db) / bg_ver / '.'.join([smp_site.lower(), 'xlsx'])
        bg_info[smp_site] = pool_bg.apply_async(get_background_info, (bg_file, )).get()
    pool_bg.close()
    pool_bg.join()

    # 获取3, 4级致病物种列表
    pc_file = Path(db) / 'PathogenicClassificationList.txt'
    df_pc = pd.read_csv(pc_file, sep='\t', dtype=str, header=1)
    df_pc = df_pc[(df_pc['致病综合打分'] == '3') | (df_pc['致病综合打分'] == '4')]
    target_pathogen = set(df_pc['l_name_s'])

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
            pool_fill.apply_async(fill, args=(tax_file, normalize, qc_file, cov_file, sub_bg_info, nc_spe2tag2rc,
                                              spe2db_code, row[1], row[2], target_pathogen))
    pool_fill.close()
    pool_fill.join()

    # 批次统计：共有物种检出
    spe2smp2comment, df_merge, spe2kingdom = stat_common(df_smp, workdir, target_pathogen)
    if len(spe2smp2comment) == 0:
        print('Warning: Common Stat Data is empty, please check it.')

    # Fill comment of common species
    for file in 'Total_Detail.txt', 'Pathogeny_Detail.txt':
        for row in df_smp.itertuples():
            tax_file = plg_dir / row[1] / file
            if not tax_file.exists():
                print(f'Error: {str(tax_file)} not exists, please check it.')
                continue
            else:
                fill_common(tax_file, row[1], spe2smp2comment)

    # 指定强阳物种统计图输出路径
    for smp in df_smp[0]:
        batch_plot_dir = plg_dir / smp / 'Batch_species_stat_plots'
        if not batch_plot_dir.exists():
            os.system(f'mkdir -p {str(batch_plot_dir)}')

    # 绘制强阳物种统计图
    pool_plot = Pool(int(thread))
    df_sub = df_merge[list(df_smp[0])].astype(int)
    x_data = list(df_sub.columns)
    for row in df_sub.itertuples():
        kd = spe2kingdom.get(row[0])
        if not kd:
            continue
        y_data = row[1:]
        if max(y_data) < 1000:
            continue
        # name = re.sub(r'[^a-zA-Z0-9]+', '_', row[0])
        name = re.sub(r'[ /]', '_', row[0])
        for smp in x_data:
            outfile = plg_dir / smp / 'Batch_species_stat_plots' / f'{name}.batch_stat.png'
            pool_plot.apply_async(draw_plot, args=(x_data, y_data, name, smp, outfile))
    pool_plot.close()
    pool_plot.join()


if __name__ == '__main__':
    main()
