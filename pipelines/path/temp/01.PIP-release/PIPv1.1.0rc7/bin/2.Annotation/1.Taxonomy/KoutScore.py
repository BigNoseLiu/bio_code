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
    """ 获取数据库中物种的注释信息 """
    tax_file = Path(db) / 'all.Taxonomy.txt'
    check_file(tax_file)

    tax_infos = defaultdict(defaultdict)

    df = pd.read_csv(tax_file, sep='\t', dtype=str)
    for row in df.itertuples():
        if row[5] == 'Viruses':
            tax_infos['taxid2name'][row[3]] = row[13]
        else:
            tax_infos['taxid2name'][row[3]] = row[11]
            tax_infos['taxid2name'][row[4]] = row[11]
        tax_infos['seqid2speid'][row[3]] = row[4]
        tax_infos['seqid2speid'][row[4]] = row[4]
        tax_infos['taxid2kingdom'][row[3]] = row[5]
        tax_infos['taxid2kingdom'][row[4]] = row[5]

    # Group infos for target groups
    grp_file = Path(db) / 'GroupList.info'
    check_file(grp_file)
    spe_tid2grp_tid = {}
    df_grp = pd.read_csv(grp_file, sep='\t', dtype=str)
    df_grp = df_grp[df_grp['Focus'] == '*']
    for row in df_grp[['GroupTaxid', 'SpeciesTaxid']].itertuples():
        spe_tid2grp_tid[row[2]] = row[1]

    # group to species set
    grp2spe_set = defaultdict(set)
    for row in df_grp[['GroupName', 'SpeciesName', 'GroupTaxid']].itertuples():
        grp2spe_set[row[1]].add(row[2])
        tax_infos['taxid2name'][row[3]] = row[1]
    tax_infos['spe_tid2grp_tid'] = spe_tid2grp_tid

    return tax_infos


def get_fail(df, col, num, set_obj):
    """ 根据col获取病毒/非病毒需要输出进行blast验证的比对结果，每个物种取前num个 """
    df_sub = df[df[col].isin(set_obj)]
    for value in set(df_sub[col].tolist()):
        df_temp = df_sub[df_sub[col] == value]
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
            if ll[2] in tax_info['taxid2kingdom']:
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

    name2rates = defaultdict(list)
    data = []
    with open(kout, 'r') as rfp:
        for _, line in enumerate(rfp):
            if line.startswith('U'):
                continue
            ll = line.strip().split('\t')
            if ll[2] == '9606' or ll[4].startswith('9606:') or ' 9606:' in ll[4]:
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
            kingdom = tax_info['taxid2kingdom'].get(s_tid)
            if kingdom == 'Viruses':
                if s_tid == ll[2]:
                    # 统计病毒种水平
                    for mapping in ll[4].split(' '):
                        m = mapping.split(':')
                        total_kmer += int(m[1])
                        if tax_info['seqid2speid'].get(m[0]) in tax_info['taxid2name']:
                            name = tax_info['taxid2name'].get(tax_info['seqid2speid'].get(m[0]))
                            spe_kmer[name] += int(m[1])
                else:
                    for mapping in ll[4].split(' '):
                        m = mapping.split(':')
                        total_kmer += int(m[1])
                        # 病毒亚种的kmer数计算包括该亚种及该亚种对应的种水平
                        if m[0] == s_tid:
                            name = tax_info['taxid2name'].get(ll[2])
                            spe_kmer[name] += int(m[1])
                        elif m[0] in tax_info['taxid2name']:
                            name = tax_info['taxid2name'].get(m[0])
                            spe_kmer[name] += int(m[1])
            else:
                if grp_tid:
                    for mapping in ll[4].split(' '):
                        m = mapping.split(':')
                        total_kmer += int(m[1])
                        tmp_stid = tax_info['seqid2speid'].get(m[0])
                        if tmp_stid in tax_info['spe_tid2grp_tid']:
                            name = tax_info['taxid2name'].get(tax_info['spe_tid2grp_tid'].get(tmp_stid))
                            spe_kmer[name] += int(m[1])
                        elif m[0] == grp_tid:
                            name = tax_info['taxid2name'].get(grp_tid)
                            spe_kmer[name] += int(m[1])
                else:
                    for mapping in ll[4].split(' '):
                        m = mapping.split(':')
                        total_kmer += int(m[1])
                        if m[0] in tax_info['taxid2name']:
                            name = tax_info['taxid2name'].get(m[0])
                            spe_kmer[name] += int(m[1])
                        # elif m[0] == s_tid:
                        #     name = tax_info['taxid2name'].get(s_tid)
                        #     spe_kmer[name] += int(m[1])
            # filtered reads that are multi mapping
            if len(spe_kmer) != 1:
                continue
            rate = list(spe_kmer.values())[0] / total_kmer
            data.append(ll + [list(spe_kmer.values())[0], rate])
            name = list(spe_kmer.keys())[0]
            name2rates[name].append(rate)

    # 计算每个物种对应的最高kmer score
    # 高于cutoff输出最高kmer score, 低于cutoff则输出比对结果
    score_out = Path(out).parent / 'Species-HighConfidenceKmerScore.txt'
    score_out_fp = score_out.open('w')
    spe_fail = set()
    for name, rates in name2rates.items():
        if not rates:
            continue
        if max(rates) >= float(kmercut):
            score_out_fp.write('\t'.join([name, "{:.2f}".format(max(rates) * 100)]) + '\n')
        else:
            spe_fail.add(name)
    score_out_fp.close()

    # 输出需要进行blast验证的序列比对结果
    with open(out, 'w') as wfp:
        df = pd.DataFrame(data, columns=None)
        df['kingdom'] = [tax_info['spe2kingdom'].get(i) for i in df[2]]
        name_lst = []
        for row in df.itertuples():
            if row[8] == 'Viruses':
                name = tax_info['taxid2name'].get(row[3])
            else:
                name = tax_info['taxid2name'].get(row[3])
            name_lst.append(name)
        df['name'] = name_lst

        for row in get_fail(df, 'name', int(top), spe_fail):
            wfp.write(row)
        #
        # # get alignment for validating from viral
        # df_viral = df[df['kingdom'] == 'Viruses']
        # if not df_viral.empty:
        #     for row in get_fail(df_viral, 's_tid', int(top), tax_info['seqid2speid'], spe_fail):
        #         wfp.write(row)
        #     for row in get_fail(df_viral, 's1_tid', int(top), tax_info['seqid2s1_tid'], spe_fail):
        #         wfp.write(row)
        # # get alignment for validating from non viral
        # df_non_viral = df[df['kingdom'] != 'Viruses']
        # if not df_non_viral.empty:
        #     for row in get_fail(df_non_viral, 's_tid', int(top), tax_info['seqid2speid'], spe_fail):
        #         wfp.write(row)
        #     for row in get_fail(df_non_viral, 'grp_tid', int(top), tax_info['spe_tid2grp_tid'], spe_fail):
        #         wfp.write(row)


if __name__ == '__main__':
    main()
