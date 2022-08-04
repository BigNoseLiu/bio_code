#!/usr/bin/env python
import os
import sys

import click
import warnings
import pandas as pd
from pathlib import Path
from collections import defaultdict

warnings.filterwarnings('ignore')

############################ developer info ##############################
__AUTHOR__ = "Sujiawei"
__EMAIL__ = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__ = "2022-04-25"
##########################################################################


def zero():
    return 0


def check_file(path):
    if not Path(path).exists():
        raise f"{str(path)} not exists!"


def get_tax_info(db):
    """ 获取物种taxid相关注释信息 """
    infile1 = db / 'all.Taxonomy.txt'
    infile2 = db / 'seqTaxid2speciesTaxid'
    infile3 = db / 'seqTaxid2S1LevelTaxid'
    infile4 = db / 'GroupList.info'
    for file in infile1, infile2, infile3:
        check_file(file)

    df1 = pd.read_csv(infile1, sep='\t', dtype=str)
    df_sub = df1[df1['Focus'] == '*']
    spe2kingdom = {}
    for row in df_sub.itertuples():
        spe2kingdom[row[1]] = row[3]
        spe2kingdom[row[2]] = row[3]

    seqid2speid = {}
    with open(infile2, 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.strip().split('\t')
            if ll:
                if ll[1] != 'NULL':
                    seqid2speid[ll[0]] = ll[1]
                    if ll[0] not in spe2kingdom and ll[1] in spe2kingdom:
                        spe2kingdom[ll[0]] = spe2kingdom.get(ll[1])

    seqid2s1_tid = {}
    with open(infile3, 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.strip().split('\t')
            if ll:
                seqid2s1_tid[ll[0]] = ll[1]
                seqid2speid[ll[0]] = ll[2]
                if ll[0] not in spe2kingdom and ll[2] in spe2kingdom:
                    spe2kingdom[ll[0]] = spe2kingdom.get(ll[2])

    spe_tid2grp_tid = {}
    df_grp = pd.read_csv(infile4, sep='\t', dtype=str)
    df_grp = df_grp[df_grp['Focus'] == '*']
    for row in df_grp[['GroupTaxid', 'SpeciesTaxid']].itertuples():
        spe_tid2grp_tid[row[2]] = row[1]

    return {
        'spe2kingdom': spe2kingdom,
        'seqid2speid': seqid2speid,
        'seqid2s1_tid': seqid2s1_tid,
        'spe_tid2grp_tid': spe_tid2grp_tid
    }


def get_fail(df, col, num, dict_obj, set_obj):
    """ 根据col获取病毒/非病毒需要输出进行blast验证的比对结果，每个物种取前num个 """
    if col != 'grp_tid':
        df[col] = df.apply(lambda x: dict_obj.get(x[2]), axis=1)
    else:
        # 复合群水平
        df = df[df[2].isin(set(dict_obj.values()))]
        df[col] = df[2]

    df = df[df[col].isin(set_obj)]
    for tid in set(df[col].tolist()):
        df_temp = df[df[col] == tid]
        # 按kmer score由大到小排序
        df_temp = df_temp.sort_values(by=6, ascending=False)
        for row in df_temp.iloc[:num].itertuples():
            yield '\t'.join(row[1:6]) + '\n'


@click.command(no_args_is_help=True)
@click.option('--kout', help="Kraken2 output file")
@click.option('--top', help="Mark top n reads to validate [5]", default=5)
@click.option('--db', help="Database version of PIP.")
@click.option('--out', help="Output file, tsv format")
@click.option('--kmerCut', help="kmer cut off. [0.6875]", default=0.6875)
def main(kout, top, db, out, kmercut):
    """"""
    # 检查输入文件是否存在
    check_file(kout)

    # 获取物种注释信息
    tax_info = get_tax_info(Path(db))
    target_groups = set(tax_info['spe_tid2grp_tid'].values())

    # 若没有reads比对到微生物则生成空结果
    with open(kout, 'r') as rfp:
        mapped_reads = 0
        for _, line in enumerate(rfp):
            if not line.startswith('C'):
                continue
            ll = line.split('\t')
            if ll[2] in tax_info['seqid2speid']:
                mapped_reads += 1
                break
            elif ll[2] in target_groups:
                mapped_reads += 1
                break
        if mapped_reads == 0:
            os.system(f'touch {out}')
            score_out = Path(out).parent / 'Species-HighConfidenceKmerScore.txt'
            os.system(f'touch {str(score_out)}')
            return

    tax_id2rates = defaultdict(list)
    data = []
    with open(kout, 'r') as rfp:
        for _, line in enumerate(rfp):
            if line.startswith('U'):
                continue
            ll = line.strip().split('\t')
            if ll[2] == '9606':
                continue
            s_tid = tax_info['seqid2speid'].get(ll[2])
            if ll[2] in target_groups:
                grp_tid = ll[2]
            else:
                grp_tid = False
            # 只统计种及种以下分类水平的kmer score, 重点复合群除外
            if not s_tid and not grp_tid:
                continue
            total_kmer = 0          # 单条read的kmer总数
            spe_kmer = defaultdict(zero)            # 鉴定到目标物种的kmer数
            kingdom = tax_info['spe2kingdom'].get(s_tid)
            s1_tid = tax_info['seqid2s1_tid'].get(ll[2])
            if kingdom == 'Viruses' and s1_tid:
                for mapping in ll[4].split(' '):
                    m = mapping.split(':')
                    total_kmer += int(m[1])
                    # 病毒亚种的kmer数计算包括该亚种及该亚种对应的种水平
                    if m[0] in tax_info['seqid2s1_tid']:
                        spe_kmer[tax_info['seqid2s1_tid'].get(m[0])] += int(m[1])
                    elif m[0] == s_tid:
                        spe_kmer[s1_tid] += int(m[1])
            else:
                if grp_tid:
                    for mapping in ll[4].split(' '):
                        m = mapping.split(':')
                        total_kmer += int(m[1])
                        tmp_stid = tax_info['seqid2speid'].get(m[0])
                        if tmp_stid in tax_info['spe_tid2grp_tid']:
                            spe_kmer[tax_info['spe_tid2grp_tid'].get(tmp_stid)] += int(m[1])
                        elif m[0] == grp_tid:
                            spe_kmer[grp_tid] += int(m[1])
                else:
                    for mapping in ll[4].split(' '):
                        m = mapping.split(':')
                        total_kmer += int(m[1])
                        if m[0] in tax_info['seqid2speid']:
                            spe_kmer[tax_info['seqid2speid'].get(m[0])] += int(m[1])
                        elif m[0] == s_tid:
                            spe_kmer[s_tid] += int(m[1])
            # filtered reads that are multi mapping
            if len(spe_kmer) != 1:
                continue
            rate = list(spe_kmer.values())[0] / total_kmer
            data.append(ll + [list(spe_kmer.values())[0], rate])
            if kingdom == 'Viruses' and s1_tid:
                tax_id2rates[s1_tid].append(rate)
            elif grp_tid:
                tax_id2rates[grp_tid].append(rate)
            else:
                tax_id2rates[s_tid].append(rate)

    # 计算每个物种对应的最高kmer score
    # 高于cutoff输出最高kmer score, 低于cutoff则输出比对结果
    score_out = Path(out).parent / 'Species-HighConfidenceKmerScore.txt'
    score_out_fp = score_out.open('w')
    spe_fail = set()
    for tax_id, rates in tax_id2rates.items():
        if max(rates) >= float(kmercut):
            score_out_fp.write('\t'.join([tax_id, "{:.2f}".format(max(rates) * 100)]) + '\n')
        else:
            spe_fail.add(tax_id)
    score_out_fp.close()

    # 输出需要进行blast验证的序列比对结果
    with open(out, 'w') as wfp:
        df = pd.DataFrame(data, columns=None)
        df['kingdom'] = df.apply(lambda x: tax_info['spe2kingdom'].get(x[2]), axis=1)
        # get alignment for validating from viral
        df_viral = df[df['kingdom'] == 'Viruses']
        if not df_viral.empty:
            for row in get_fail(df_viral, 's_tid', int(top), tax_info['seqid2speid'], spe_fail):
                wfp.write(row)
            for row in get_fail(df_viral, 's1_tid', int(top), tax_info['seqid2s1_tid'], spe_fail):
                wfp.write(row)
        # get alignment for validating from non viral
        df_non_viral = df[df['kingdom'] != 'Viruses']
        if not df_non_viral.empty:
            for row in get_fail(df_non_viral, 's_tid', int(top), tax_info['seqid2speid'], spe_fail):
                wfp.write(row)
            for row in get_fail(df_non_viral, 'grp_tid', int(top), tax_info['spe_tid2grp_tid'], spe_fail):
                wfp.write(row)


if __name__ == '__main__':
    main()
