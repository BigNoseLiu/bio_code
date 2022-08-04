#!/usr/bin/env python
import os
import re
import sys
import json
import math
import time

import click
import warnings
import subprocess
import pandas as pd
import multiprocessing
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as pat
from multiprocessing import Pool
from collections import defaultdict
Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg
from configure import LogExceptions

warnings.filterwarnings('ignore')

############################ developer info ##############################
__AUTHOR__ = "Sujiawei"
__EMAIL__ = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__ = "2022-07-15"
##########################################################################

KINGDOM2TAG = {
    'Bacteria': 'Bacteria',
    'Eukaryota:Fungi': 'Eukaryota_Fungi',
    'Eukaryota:Parasite': 'Eukaryota_Parasite',
    'Eukaryota:Protozoa': 'Eukaryota_Protozoa',
    'SpecialPathogens': 'SpecialPathogens',
    'Viruses': 'Viruses'
}


def check_file(path):
    if not Path(path).exists():
        print(f'{str(path)} not exists! Please check it again.')


def create_dir(path_list):
    for path in path_list:
        if not Path(path).exists():
            os.system(f'mkdir -p {str(path)}')


############################### Annotating by DB ###############################
def get_tax_info(db):
    """ 获取数据库中物种的注释信息 """
    tax_file = Path(db) / 'all.Taxonomy.txt'
    check_file(tax_file)

    name2tid = {}
    tid2spe_tid = {}
    # tid2spe_len = {}
    name2spe_len = {}
    spe_tid2grp_tid = {}
    grp2spe_tid_set = defaultdict(set)

    df = pd.read_csv(tax_file, sep='\t', dtype=str)
    for row in df.itertuples():
        name2tid[row[10]] = row[17]
        if row[5] == 'Viruses':
            name2tid[row[13]] = row[1]
            tid2spe_tid[row[3]] = row[1]
            # tid2spe_len[row[3]] = int(row[14])
            name2spe_len[row[13]] = int(row[14])
            if row[3] == row[4]:
                tid2spe_tid[row[4]] = row[1]
                # tid2spe_len[row[4]] = int(row[15])
                name2spe_len[row[13]] = int(row[15])
        else:
            name2tid[row[11]] = row[4]
            tid2spe_tid[row[3]] = row[4]
            tid2spe_tid[row[4]] = row[4]
            # tid2spe_len[row[4]] = int(row[15])
            name2spe_len[row[11]] = int(row[15])

    # Group infos for target groups
    grp_file = Path(db) / 'GroupList.info'
    check_file(grp_file)
    df_grp = pd.read_csv(grp_file, sep='\t', dtype=str)
    df_grp = df_grp[df_grp['Focus'] == '*']
    for row in df_grp.itertuples():
        name2tid[row[2]] = row[1]
        grp2spe_tid_set[row[2]].add(row[3])
        spe_tid2grp_tid[row[3]] = row[1]
        if row[2] not in name2spe_len:
            # tid2spe_len[row[1]] = tid2spe_len.get(row[3], 0)
            name2spe_len[row[2]] = name2spe_len.get(row[4], 0)
        else:
            # tid2spe_len[row[1]] += tid2spe_len.get(row[3], 0)
            name2spe_len[row[2]] += name2spe_len.get(row[4], 0)

    return {
        'name2tid': name2tid,
        'tid2spe_tid': tid2spe_tid,
        'grp2spe_tid_set': grp2spe_tid_set,
        'grp_set': set(df_grp['GroupName']),
        'grp_tid_set': set(df_grp['GroupTaxid']),
        'spe_tid2grp_tid': spe_tid2grp_tid,
        # 'tid2spe_len': tid2spe_len,
        'name2spe_len': name2spe_len,
    }


def align(cov_db, target_tid2old_kd, kingdom_tag, out_dir, thread):
    """ kraken2 比对物种参考基因组序列 """
    fna_lst = []
    for tid, old_kd in target_tid2old_kd.items():
        spe_fna = Path(cfg.coverage) / 'v1.0.1' / old_kd / f'{tid}.fna'
        if spe_fna.exists():
            fna_lst.append(str(spe_fna))
    query_seq = ' '.join(fna_lst)
    krk_out = out_dir / f'{kingdom_tag}.kout'
    cmd = f'{cfg.kraken} --db {str(cov_db)} --threads {thread} <(cat {query_seq}) > {krk_out}'
    subprocess.call(['bash', '-c', cmd])


def draw_plot(data, name, outfile):
    """ draw coverage plot """
    fig, ax = plt.subplots()

    x_data = [int(i) for i in data['Unique'].keys()]
    y1_data = list(data['Unique'].values())
    y2_data = list(data['Multiple'].values())
    y3_data = list(data['AvgDepth'].values())

    # 左轴
    ax.bar(x_data, y2_data, color='darkgrey')
    ax.bar(x_data, y1_data, color='dodgerblue')
    if data['GenomeLength'] < 10 ** 4:
        ax.set_xlabel('Position on the genome (bp)', fontdict={'fontsize': 14})
    else:
        ax.set_xlabel('Position on the genome (Mbp)', fontdict={'fontsize': 14})
    ax.set_ylabel('Mapping Reads', fontdict={'fontsize': 14})
    y1_ticks = [0] + [math.ceil(i * max(y1_data + y2_data) / 5) for i in range(1, 6, 1)]
    ax.set_yticks(y1_ticks)
    ax.tick_params(axis='y', colors='dodgerblue')

    # ax.fill_between(x_data, y2_data, facecolor='darkgrey')
    # ax.fill_between(x_data, y1_data, facecolor='dodgerblue')

    # 处理边距
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('dodgerblue')
    ax.spines['bottom'].set_position(('data', 0))
    ax.margins(False)

    # 右轴
    ax2 = ax.twinx()
    ax2.plot(x_data, y3_data, color='darkorange', label='Average Depth')
    ax2.set_xlabel('Position on the genome (Mbp)', fontdict={'fontsize': 14})
    ax2.set_ylabel('Average Depth', fontdict={'fontsize': 14})
    y3_max = max(y3_data)
    ax2.set_ylim(min(y3_data), y3_max)
    y2_sep = y3_max / 5
    y2_ticks = [0] + [round(i * y2_sep, 4) for i in range(1, 5, 1)] + [round(y3_max, 4)]
    ax2.set_yticks(y2_ticks)
    ax2.tick_params(axis='y', colors='darkorange')

    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_color('darkorange')
    ax2.margins(False)

    # 自定义横坐标轴刻度标签
    sep = len(x_data) / 10
    bin_len = data['GenomeLength'] / 10
    x_idx = [0] + [i for i in x_data if i % int(sep) == 0]
    x_ticks = [0]
    if data['GenomeLength'] < 10 ** 4:
        for i in range(1, 11, 1):
            if i % 2 == 1:
                x_ticks.append('')
            else:
                x_ticks.append(str(math.ceil(i * bin_len)))
    else:
        for i in range(1, 11, 1):
            if i % 2 == 1:
                x_ticks.append('')
            else:
                x_ticks.append('{:.2f}'.format(i * bin_len / 10 ** 6))
    plt.xticks(x_idx[:11], x_ticks)

    # 处理图例、标题
    plt.subplots_adjust(top=0.77, bottom=0.14)
    if len(name) > 45:
        fontsize = 12
    elif len(name) > 35:
        fontsize = 17
    else:
        fontsize = 19
    plt.suptitle(f"{name}", fontsize=fontsize, fontweight='semibold', y=0.95)
    # plt.xlabel('Position on the genome (Mbp)', fontdict={'fontsize': 14})
    # plt.ylabel('Mapping Reads', fontdict={'fontsize': 14})
    patches = [pat.Patch(color='dodgerblue', label='Unique Mapping'),
               pat.Patch(color='darkgrey', label='Multiple Mapping'),
               pat.Patch(color='darkorange', label='Average Depth')]
    plt.gca().legend(handles=patches, loc='upper center', bbox_to_anchor=(0.55, 1.19), ncol=3, fontsize=9)

    # 调整图片位置
    plt.gcf().subplots_adjust(right=0.829)

    plt.savefig(outfile)


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


def stat4cov(species, spe_reads, genus, gns_reads, read_len, tax_info, acc_info, df_data, out_png):
    """ 根据kraken2输出的out文件统计比对结果 """
    spe_tid = tax_info['name2tid'].get(species, '0')
    gns_tid = tax_info['name2tid'].get(genus, '0')
    if gns_tid:
        gns_tid = str(gns_tid).split('.')[0]

    df_kout = df_data[df_data['SpeTaxid'] == spe_tid]
    if df_kout.empty:
        print('################')
        print(f'\"{species}\" mapped nothing. SpeciesTaxid={species}, SpeciesReads={spe_reads}')

        return

    seq_id2accession = acc_info.get(spe_tid)
    if seq_id2accession:
        asm_lst = []
        for seqid in df_kout[1]:
            asm = seq_id2accession.get(seqid)
            if asm and '-' not in asm:
                asm_lst.append(asm)
            else:
                asm_lst.append(seqid)
        df_kout['Assembly'] = asm_lst
    else:
        df_kout['Assembly'] = df_kout[1]

    # 统计覆盖区总长
    k_cnt_lst = []
    if species in tax_info['grp_set']:
        for alignment in df_kout[4]:
            k_cnt = 0
            for m in alignment.split(' '):
                m = m.split(':')
                if m in ['0', '1', '2', 'A']:
                    continue
                tid = tax_info['tid2spe_tid'].get(m[0])
                if tax_info['spe_tid2grp_tid'].get(tid) == spe_tid:
                    k_cnt += int(m[1])
            k_cnt_lst.append(k_cnt)
    else:
        for alignment in df_kout[4]:
            k_cnt = 0
            for m in alignment.split(' '):
                m = m.split(':')
                if m in ['0', '1', '2', 'A']:
                    continue
                if tax_info['tid2spe_tid'].get(m[0]) == spe_tid:
                    k_cnt += int(m[1])
            k_cnt_lst.append(k_cnt)
    df_kout['k_cnt'] = k_cnt_lst

    # 获取最优比对参考基因组
    if len(set(df_kout['Assembly'])) != 1:
        df_gb = df_kout.groupby('Assembly').sum()
        max_k_cnt = max(df_gb['k_cnt'])
        for row in df_gb.itertuples():
            if row[1] == max_k_cnt:
                df_kout = df_kout[df_kout['Assembly'] == row[0]]
                break

    mapped_pos_len = sum(df_kout['k_cnt'])
    mapped_gnm_len = sum(df_kout[3].astype(int))
    # genome_len = tax_info['tid2spe_len'].get(spe_tid, mapped_pos_len)
    genome_len = tax_info['name2spe_len'].get(species, mapped_pos_len)
    bin_num = 100
    bin_len = math.ceil(mapped_gnm_len / bin_num)

    result = {
        # 'Coverage': '{:.3f}'.format((mapped_pos_len + 34) / max(mapped_gnm_len, genome_len) * 100),
        'Unique': {i + 1: 0 for i in range(bin_num)},
        'Multiple': {i + 1: 0 for i in range(bin_num)},
        # 'AvgDepth': {i + 1: 0 for i in range(bin_num)},
    }

    pos = 0
    u_idx2len = {i + 1: 0 for i in range(bin_num)}     # unique mapped reads
    m_idx2len = {i + 1: 0 for i in range(bin_num)}     # multiple mapped reads
    if species in tax_info['grp_set']:
        # 复合群
        for alignment in df_kout[4]:
            for m in alignment.strip().split(' '):
                m = m.split(':')
                kmer = int(m[1])
                pos += kmer
                tid = tax_info['tid2spe_tid'].get(m[0])
                if tax_info['spe_tid2grp_tid'].get(tid) == spe_tid:
                    u_idx2len[math.ceil(pos / bin_len)] += kmer
                if m[0] == gns_tid:
                    m_idx2len[math.ceil(pos / bin_len)] += kmer
    else:
        for alignment in df_kout[4]:
            for m in alignment.strip().split(' '):
                m = m.split(':')
                kmer = int(m[1])
                pos += kmer
                if tax_info['tid2spe_tid'].get(m[0]) == spe_tid:
                    u_idx2len[math.ceil(pos / bin_len)] += kmer
                if m[0] == gns_tid:
                    m_idx2len[math.ceil(pos / bin_len)] += kmer

    if not (mapped_pos_len * spe_reads):
        print(f'Warning: species={species}, mapped_pos_len={mapped_pos_len}, spe_reads={spe_reads}, ZeroDivisionError')
        return {
            'Coverage': round(get_quotient(spe_reads * read_len, genome_len), 3),
            'BinCoverage': '-',
            'GenomeLength': genome_len,
            'Unique': {},
            'Multiple': {},
            'AvgDepth': {},
        }

    for k, v in u_idx2len.items():
        interval_rc = int(v / mapped_pos_len * spe_reads)
        result['Unique'][k] = interval_rc
    if sum(result['Unique'].values()) < spe_reads:
        diff_value = spe_reads - sum(result['Unique'].values())
        for k, v in u_idx2len.items():
            if v > 0:
                result['Unique'][k] += 1
                diff_value -= 1
                if diff_value == 0:
                    break

    multi_mapping_len = sum(m_idx2len.values())

    if gns_reads and multi_mapping_len:
        for k, v in m_idx2len.items():
            interval_rc = int(v / multi_mapping_len * gns_reads)
            result['Multiple'][k] = interval_rc
        if sum(result['Multiple'].values()) < gns_reads:
            diff_value = gns_reads - sum(result['Multiple'].values())
            for k, v in m_idx2len.items():
                if v > 0:
                    result['Multiple'][k] += 1
                    diff_value -= 1
                    if diff_value == 0:
                        break

    # 重定义输出数据
    ret_result = {
        'Coverage': '{:.3f}'.format((mapped_pos_len + 34) / max(mapped_gnm_len, genome_len) * 100),
        'BinCoverage': math.floor(len([i for i in result['Unique'].values() if i]) / bin_num * 100),
        'GenomeLength': genome_len,
        'Unique': {},
        'Multiple': {},
        'AvgDepth': {},
    }
    new_bin_len = math.ceil(mapped_gnm_len / 50)
    for i in range(1, 101, 2):
        new_idx = int((i + 1) / 2)
        new_uni = result['Unique'][i] + result['Unique'][i+1]
        new_mul = result['Multiple'][i] + result['Multiple'][i+1]
        ret_result['Unique'][new_idx] = new_uni
        ret_result['Multiple'][new_idx] = new_mul
        # 统计Bin上的平均覆盖深度
        ret_result['AvgDepth'][new_idx] = round((new_uni + new_mul) * read_len / new_bin_len, 4)

    draw_plot(ret_result, species, out_png)

    return ret_result


def stat(indir, tag, df_data, gns2sum_src, gns2reads, tax_info, read_len):
    """"""
    krk_out = Path(indir) / f'Coverage/02.kraken/{tag}.kout'
    if not krk_out.exists():
        print(f'Error: {krk_out} not exists, maybe this step failed at \"Build Kraken Database\"')
        return {}

    df_kout = pd.read_csv(krk_out, sep='\t', header=None, dtype=str)
    if df_kout.empty:
        return {}

    df_kout['SpeTaxid'] = [tax_info['tid2spe_tid'].get(i.split('|')[1], i.split('|')[1]) for i in df_kout[1]]

    spe2res = {}

    acc_file = Path(cfg.coverage) / 'v1.0.1' / tag / f'{tag}.acc.info.json'
    with open(acc_file, 'r') as rfp:
        acc_info = json.load(rfp)
    for row in df_data[['SpeciesSN', 'SpeciesReads', 'GenusSN']].itertuples():
        name = re.sub(r'[/?#$%^&*@!~|\\()\[\] \'\"]', '_', row[1])
        out_png = Path(indir) / f'Coverage/04.plot/{tag}' / (name + '.coverage.png')
        create_dir([out_png.parent])
        if str(row[3]) in ['', 'nan', '-']:
            genus_reads = 0
        else:
            genus_reads = int(int(row[2]) / gns2sum_src.get(row[3]) * gns2reads.get(row[3], 0))
        spe2res[row[1]] = stat4cov(row[1], int(row[2]), row[3], genus_reads, read_len,
                                   tax_info, acc_info, df_kout, out_png)

    return spe2res


@click.command(no_args_is_help=True)
@click.option('-i', '--indir', help='Input dir <str>')
@click.option('-o', '--outdir', help='Output dir <str>')
@click.option('--db', help='PIDB path <str>')
@click.option('-t', '--thread', help='How many threads can be used')
def main(indir, db, thread, outdir):
    """ 覆盖度模块：比对、统计、画图 """
    s_time = time.time()

    file_total = Path(indir) / 'Total_Detail.xlsx'
    df = pd.read_excel(file_total, dtype=str)
    if df.empty:
        print(f'Pathogen result file is empty, please check it: {str(file_total)}')
        return

    align_dir = Path(outdir) / 'Coverage/02.kraken'
    stat_dir = Path(outdir) / 'Coverage/03.stat'
    plot_dir = Path(outdir) / 'Coverage/04.plot'
    create_dir([align_dir, stat_dir, plot_dir])

    # get annotation from PIDB
    tax_info = get_tax_info(db)

    # print err log from subprocess
    multiprocessing.log_to_stderr()

    # align genomes
    for kd, tag in KINGDOM2TAG.items():
        df_sub = df[df['Kingdom'] == kd]
        if df_sub.empty:
            continue
        target_tid2old_kd = {}
        for row in df_sub[['SpeciesSN', 'rawKingdom']].itertuples():
            if tax_info['name2tid'].get(row[1]):
                target_tid2old_kd[tax_info['name2tid'].get(row[1])] = KINGDOM2TAG.get(row[2])
        cov_db = Path(indir) / 'Coverage/01.db_path'
        align(cov_db, target_tid2old_kd, tag, align_dir, str(thread))

    with open(Path(outdir) / 'Coverage/01.db_path/library/library.fna', 'r') as r:
        next(r)
        read_len = len(next(r).strip())

    # 获取实际比对到genus水平的reads数
    df_temp = df[['GenusSN', 'SpeciesReads']].astype({'SpeciesReads': int})
    df_gb = df_temp.groupby('GenusSN').sum()
    # genus: sum(species reads count)
    gns2sum_src = {row[0]: row[1] for row in df_gb.itertuples()}
    gns2reads = {}
    for row in df[['GenusSN', 'GenusReads']].itertuples():
        gns2reads[row[1]] = int(row[2]) - gns2sum_src.get(row[1], 0)
        # some viruses may have no Genus
        if gns2reads[row[1]] < 0:
            gns2reads[row[1]] = 0

    # limit pool threads num
    stat_pool = Pool(int(thread) - 2)
    idx2res = {}
    idx = 0
    for kd, tag in KINGDOM2TAG.items():
        df_sub = df[df['Kingdom'] == kd]
        if df_sub.empty:
            continue
        step_num = 50
        if len(df_sub) > step_num:
            for i in range(int(len(df_sub) / step_num) + 1):
                df_tmp = df_sub.iloc[i * step_num:(i + 1) * step_num, :]
                idx2res[idx] = stat_pool.apply_async(LogExceptions(stat),
                                                args=(indir, tag, df_tmp, gns2sum_src, gns2reads, tax_info, read_len))
                idx += 1
        else:
            idx2res[idx] = stat_pool.apply_async(LogExceptions(stat),
                                                 args=(indir, tag, df_sub, gns2sum_src, gns2reads, tax_info, read_len))
            idx += 1
    stat_pool.close()
    stat_pool.join()

    spe2stat = {}
    for idx, res in idx2res.items():
        spe2stat.update(res.get())

    with open(stat_dir / 'stat.json', 'w') as wfp:
        json.dump(spe2stat, wfp, indent=4)

    print(f'Coverage Elapsed Time (s): ', (time.time() - s_time) / 60, 'minutes.')


if __name__ == '__main__':
    main()
