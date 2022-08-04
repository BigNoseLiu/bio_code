#!/usr/bin/env python
import os
import re
import sys
import json
import math
import click
import warnings
import traceback
import subprocess
import pandas as pd
from pathlib import Path
import multiprocessing
from multiprocessing import Pool
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as pat

Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../conf'))
import configure as cfg
from configure import LogExceptions

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
    try:
        numerator = int(numerator)
        denominator = int(denominator)
    except:
        print(f'{numerator}, {denominator} can not convert to int, please check it.')
        return 0
    if denominator == 0:
        return 0
    else:
        return numerator / denominator


def stat_common(df_sample, workdir, target_pathogen):
    """ 统计共有检出 """
    # 合并批次所有样本
    df_merge = pd.DataFrame()  # 未均一化结果
    df_merge_normal = pd.DataFrame()  # 已均一化结果
    spe2kingdom = {}
    spe2species_cn = {}
    for row in df_sample.itertuples():
        infile = Path(workdir) / 'plugin_out/Result' / row[1] / 'Total_Detail.raw.txt'
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
    spe2smp2label = defaultdict(dict)
    spe2smp2comment = defaultdict(dict)
    spe2comment = {}
    total_smp = len(df_sample)

    for row in df_merge_normal.itertuples():
        # 统计样本物种检出与批次物种检出最大值的比值
        batch_max = max(row[1:])
        for i, v in enumerate(row[1:]):
            ratio = '({:.2f})'.format(get_quotient(v, batch_max))
            spe2smp2ratio[row[0]][idx2smp.get(i)] = ratio
            spe2smp2label[row[0]][idx2smp.get(i)] = ''.join([str(int(v)), '/', str(int(batch_max)), ' ', ratio])

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
            ratio = spe2smp2ratio[row[0]].get(smp, '-')
            spe2smp2comment[row[0]][smp] = f'{common}|{str(value2rank.get(v))}{ratio}'

    # 补充物种注释信息
    smp_cols = list(df_merge.columns)
    out_cols = ['Kingdom', 'SpeciesSN', 'SpeciesCN', 'Comment'] + smp_cols
    for df in df_merge, df_merge_normal:
        df['Kingdom'] = [spe2kingdom.get(i, '-') for i in df.index]
        df['SpeciesSN'] = df.index
        df['SpeciesCN'] = [spe2species_cn.get(i, '-') for i in df.index]
        df['Comment'] = [spe2comment.get(i, '-') for i in df.index]

    for col in (set(df_sample[0]) ^ set(smp_cols)):
        df_merge[col] = 0
        df_merge_normal[col] = 0

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

    return spe2smp2comment, df_merge, spe2kingdom, spe2smp2label


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
        # new rule
        if re.match(r'N\d+H', row[1]):
            # DNA 阴性质控样本-去宿主处理
            tag2smp_lib['D_nc_H'] = '_'.join([row[1], row[2]])
        elif re.match(r'N\d+D', row[1]):
            # DNA 阴性质控样本
            tag2smp_lib['D_nc'] = '_'.join([row[1], row[2]])
        elif re.match(r'N\d+R', row[1]):
            # RNA 阴性质控样本
            tag2smp_lib['R_nc'] = '_'.join([row[1], row[2]])
        # old rule
        elif re.match(r'Dzk.+DH', row[1]):    # old rule, DNA Negative Control, delete host
            tag2smp_lib['D_nc_H'] = '_'.join([row[1], row[2]])
        elif re.match(r'Rzk\d+R', row[1]):    # old rule, RNA Negative Control
            tag2smp_lib['R_nc'] = '_'.join([row[1], row[2]])

    if not tag2smp_lib:
        print(f"Warning: NC sample name's format is unvalidated")
        return defaultdict(dict)

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


# def sort_filter(df_data, target_pathogen, top=10):
#     """ 按界进行归类，再按属reads、种reads分别由高到低进行排序; 根据致病菌分级列表进行筛选, 每个属最多保留前top个物种结果 """
#     df_ret = pd.DataFrame(columns=df_data.columns)
#     for kd in KINGDOM2TAG.keys():
#         df_kd = df_data[df_data['Kingdom'] == kd]
#         if df_kd.empty:
#             continue
#         df_kd = df_kd.astype({'GenusReads': int, 'SpeciesReads': int})
#         df_sort = df_kd.sort_values(by=['GenusReads', 'SpeciesReads'], ascending=False)
#         df_filter = pd.DataFrame(columns=df_data.columns)
#         for gns in set(df_sort['GenusSN']):
#             df_gns = df_sort[df_sort['GenusSN'] == gns]
#             # is pathogens
#             df_is_ptg = df_gns[df_gns['SpeciesSN'].isin(target_pathogen)]
#             df_filter = pd.concat([df_filter, df_is_ptg])
#             if len(df_is_ptg) >= 10:
#                 continue
#             # is not pathogens
#             df_not_ptg = df_gns[~df_gns['SpeciesSN'].isin(target_pathogen)]
#             df_filter = pd.concat([df_filter, df_not_ptg.iloc[:top-len(df_is_ptg)]])
#         df_filter = df_filter.sort_values(by=['GenusReads', 'SpeciesReads'], ascending=False)
#         df_ret = pd.concat([df_ret, df_filter])
#     return df_ret


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
    df_data['Freq/[Min-25%-Median-75%-Max](Historically)'] = [bg_info['spe2bg_fre'].get(i, '-') for i in
                                                              df_data['SpeciesSN']]
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

    if re.match(r'Dzk.+DH', smp):     # old rule, DNA Negative Control, delete host
        tag = 'D_nc_H'
    elif re.match(r'N\d+H', smp):      # new rule, DNA Negative Control, delete host
        tag = 'D_nc_H'
    elif re.match(r'N\d+D', smp):      # new rule, DNA Negative Control
        tag = 'D_nc'
    elif re.match(r'Rzk\d+R', smp):    # old rule, RNA Negative Control
        tag = 'R_nc'
    elif re.match(r'N\d+R', smp):      # new rule, RNA Negative Control
        tag = 'R_nc'
    elif re.match(r'[A-Z]\d+H', smp):  # new rule, Clinical DNA sample, delete host
        tag = 'D_nc_H'
    elif re.match(r'[A-Z]\d+D', smp):  # new rule, Clinical DNA sample
        tag = 'D_nc'
    elif re.match(r'[A-Z]\d+R', smp):  # new rule, Clinical RNA sample
        tag = 'R_nc'
    else:
        tag = 'D_nc'

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
            nc_fdc.append('{:.3f}'.format(get_quotient(row[2], row[1])))
    df_data['FoldChange(NC)'] = nc_fdc

    return df_data


def fill_common(file_path, sample, spe2smp2comment):
    """ 补充共有检出统计 """
    df_data = pd.read_csv(file_path, sep='\t', dtype=str)

    if df_data.empty:
        print(f'Warning: {file_path} is empty, please check it.')
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
        os.system(f'mkdir -p {str(file_path.parent)}')

    # Data Input
    df_data = pd.read_csv(file_path, sep='\t', dtype=str)
    if df_data.empty:
        print(f'Warning: {file_path} is empty, please check it.')
        return

    # # Sort and filter data
    # df_data = sort_filter(df_data, target_pathogen)

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


def draw_plot(x_data, y_data, species_name, sample, spe2smp2label, outfile):
    """ 绘制批次强阳物种统计柱状图 """
    x2y = {x_data[i]: y_data[i] for i in range(len(x_data))}
    data_order = sorted(x2y.items(), key=lambda l: l[1], reverse=False)
    x = [''] + [t[0] for t in data_order]
    y = [0] + [t[1] for t in data_order]

    fig, ax = plt.subplots(figsize=(6, 14))
    color_lst = []
    smp_idx = 0
    for i, smp in enumerate(x):
        if smp == sample:
            smp_idx = i
            color_lst.append('coral')
        else:
            color_lst.append('dodgerblue')
    ax.barh(x, y, color=color_lst, )

    # 坐标轴设置
    ax.set_xlabel('Unique Mapping Reads', fontdict={'fontsize': 18})
    ax.set_ylabel('Sample ID', fontdict={'fontsize': 18})
    ax.set_xlim(0, max(y))
    ax.set_xticks([0, int(max(y) / 2), max(y)])
    plt.yticks(size=13)

    # 处理边距
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('dodgerblue')
    ax.spines['bottom'].set_position(('data', 0))
    ax.margins(False)

    # 处理图例、标题
    plt.subplots_adjust(top=0.77, bottom=0.14)
    plt.suptitle(f"{species_name}", fontsize=28, fontweight='semibold', y=0.863)
    patches = [pat.Patch(color='dodgerblue', label='Others Sample ID'),
               pat.Patch(color='coral', label='Current Sample ID')]
    plt.gca().legend(handles=patches, loc='upper center', ncol=2, fontsize=12.5, bbox_to_anchor=(0.5, 1.082))

    # 调整图片位置
    plt.gcf().subplots_adjust(wspace=80)

    # 坐标轴刻度设置
    if smp_idx:
        ax.get_yticklabels()[smp_idx].set_color('coral')

    # 柱状图值标注
    if species_name in spe2smp2label:
        smp2label = spe2smp2label[species_name]
        for i, j in zip(x, y):
            plt.text(max(y), i, '     ' + smp2label.get(i, ''), fontsize=13, verticalalignment='center', )

    plt.savefig(outfile, bbox_inches='tight', pad_inches=0.5)


def filter_by_white_list(df_data, db, spe2kingdom):
    """ 输出需要绘图的物种名单 """
    white_list_set = set()
    white_list = Path(db) / 'WhiteList.txt'
    if white_list.exists():
        with open(white_list, 'r') as rfp:
            for _, line in enumerate(rfp):
                ll = line.split('\t')
                if ll:
                    white_list_set.add(ll[0])

    spe2max_rc = {row[0]: max(row[1:]) for row in df_data.itertuples()}

    ret_set = set()
    for species in df_data.index:
        if species in white_list_set:
            ret_set.add(species)
        else:
            if spe2kingdom.get(species) == 'Bacteria' and spe2max_rc.get(species, 0) >= 1000:
                ret_set.add(species)
                continue
            if spe2kingdom.get(species) == 'Eukaryota:Fungi' and spe2max_rc.get(species, 0) > 100:
                ret_set.add(species)

    return ret_set


@click.command(no_args_is_help=True)
@click.option('-w', '--workdir', help='Where the PIP output dir is <str>')
@click.option('-s', '--samples', help='Sample list file <str>')
@click.option('--nc', help='Negative Control sample\'s name <str>', default='')
@click.option('--db', help='PIDB path <str>')
@click.option('-b', '--bg_ver', help='Background database\'s version. <str>', default='BG.v1')
@click.option('-n', '--normalize', help='n M reads to normalize. [20]', default=20)
@click.option('-t', '--thread', help='How many threads can be used. [120]', default=120)
def main(workdir, samples, nc, db, bg_ver, normalize, thread):
    """ 补充阴控、背景菌等相关数据 """
    # 样本列表
    df_smp = pd.read_csv(samples, sep='\t', dtype=str, header=None)

    # 获取阴性样本数据
    nc_spe2tag2rc = get_nc(df_smp, normalize, workdir)

    # 获取Database code
    spe2db_code = {}
    all_spe2kd = {}
    with open(Path(db) / 'all.Taxonomy.txt') as rfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                if ll[4] == 'Viruses':
                    spe2db_code[ll[12]] = ll[0]
                    all_spe2kd[ll[13]] = ll[5]
                else:
                    spe2db_code[ll[10]] = ll[0]
                    all_spe2kd[ll[11]] = ll[5]

    # print err log from subprocess
    multiprocessing.log_to_stderr()

    # 获取背景菌数据库数据
    bg_info = {}
    pool_bg = Pool(int(thread))
    for smp_site in set(df_smp[len(df_smp.columns) - 1]):
        bg_file = Path(db) / bg_ver / '.'.join([smp_site.lower(), 'xlsx'])
        bg_info[smp_site] = pool_bg.apply_async(LogExceptions(get_background_info), (bg_file,)).get()
    pool_bg.close()
    pool_bg.join()

    # 获取3, 4级致病物种列表
    pc_file = Path(db) / 'PathogenicClassificationList.txt'
    df_pc = pd.read_csv(pc_file, sep='\t', dtype=str, header=1)
    if 'l_name_s' in df_pc.columns:
        # 版本1
        df_pc = df_pc[(df_pc['致病综合打分'] == '1') | (df_pc['致病综合打分'] == '2') |
                      (df_pc['致病综合打分'] == '3') | (df_pc['致病综合打分'] == '4')]
        target_pathogen = set(df_pc['l_name_s'])
    else:
        # 版本2
        df_pc = df_pc[(df_pc['致病等级'] == '1') | (df_pc['致病等级'] == '2') |
                      (df_pc['致病等级'] == '3') | (df_pc['致病等级'] == '4')]
        target_pathogen = set(df_pc['SpeciesEN'])

    qc_dir = Path(workdir) / '01.QC/01.Fastp'
    plg_dir = Path(workdir) / 'plugin_out/Result'

    # 补充覆盖度、背景菌数据、共有检出等至结果表
    pool_fill = Pool(int(thread))
    for file in 'Total_Detail.raw.txt', 'Total_Detail.txt', 'Pathogeny_Detail.txt':
        for row in df_smp.itertuples():
            tax_file = plg_dir / row[1] / file
            if not tax_file.exists():
                print(f'Error: {tax_file} not exists, please check it.')
                continue
            qc_file = qc_dir / row[1] / row[2] / f'{row[1]}.json'
            sub_bg_info = bg_info[row[-1]]
            cov_file = Path(workdir) / f'02.Annotation/01.Taxonomy/{row[1]}/{row[2]}/Coverage/03.stat/stat.json'
            pool_fill.apply_async(LogExceptions(fill),
                                  args=(tax_file, normalize, qc_file, cov_file, sub_bg_info, nc_spe2tag2rc,
                                        spe2db_code, row[1], row[2], target_pathogen))
    pool_fill.close()
    pool_fill.join()

    # 批次统计：共有物种检出
    spe2smp2comment, df_merge, spe2kingdom, spe2smp2label = stat_common(df_smp, workdir, target_pathogen)
    if len(spe2smp2comment) == 0:
        print('Warning: Common Stat Data is empty, please check it.')

    # Fill comment of common species
    for file in 'Total_Detail.raw.txt', 'Total_Detail.txt', 'Pathogeny_Detail.txt':
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
    # 处理结果为空的样本
    for col in (set(df_smp[0]) ^ set(df_merge.columns)):
        df_merge[col] = 0
    df_sub = df_merge[list(df_smp[0])].astype(int)
    x_data = list(df_sub.columns)
    spe_white_list = filter_by_white_list(df_sub, db, all_spe2kd)
    for row in df_sub.itertuples():
        if row[0] not in spe_white_list:
            continue
        y_data = row[1:]
        # name = re.sub(r'[^a-zA-Z0-9]+', '_', row[0])
        name = re.sub(r'[/?#$%^&*@!~|\\()\[\] \'\"]', '_', row[0])
        for smp in x_data:
            batch_plot_dir = plg_dir / smp / 'Batch_species_stat_plots'
            outfile = batch_plot_dir / f'{name}.batch_stat.png'
            pool_plot.apply_async(LogExceptions(draw_plot),
                                  args=(x_data, y_data, row[0], smp, spe2smp2label, outfile))
    pool_plot.close()
    pool_plot.join()


if __name__ == '__main__':
    main()
