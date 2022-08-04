#!/usr/bin/env python
import os
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


def check_file(path):
    if not Path(path).exists():
        raise f"{str(path)} not exists!"


def get_tax_info(db):
    """ 获取物种taxid相关注释信息 """
    infile1 = db / 'all.Taxonomy.txt'
    infile2 = db / 'seqTaxid2speciesTaxid'
    infile3 = db / 'seqTaxid2S1LevelTaxid'
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

    return spe2kingdom, seqid2speid, seqid2s1_tid


def get_fail(df, col, num, dict_obj, set_obj):
    """ 根据col获取病毒/非病毒需要输出进行blast验证的比对结果，每个物种取前num个 """
    df[col] = df.apply(lambda x: dict_obj.get(x[2]), axis=1)
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
    # backup kout and shuf -hy start
    #cmd = f'mv {kout} {kout}_backup; shuf {kout}_backup -o {kout}'
    #os.system(cmd)	
    # backup kout and shuf -hy end
    
    spe2kingdom, seqid2speid, seqid2s1_tid = get_tax_info(Path(db))

    check_file(kout)

    tax_id2rates = defaultdict(list)
    data = []
    with open(kout, 'r') as rfp:
        for _, line in enumerate(rfp):
            if line.startswith('U'):
                continue
            ll = line.strip().split('\t')
            if ll[2] == '9606':
                continue
            s_tid = seqid2speid.get(ll[2])
            # 只统计种及种以下分类水平的kmer score
            if not s_tid:
                continue
            total_kmer = 0          # 单条read的kmer总数
            spe_kmer = 0            # 鉴定到目标物种的kmer数
            kingdom = spe2kingdom.get(s_tid)
            s1_tid = seqid2s1_tid.get(ll[2])
            if kingdom == 'Viruses' and s1_tid:
                for mapping in ll[4].split(' '):
                    m = mapping.split(':')
                    total_kmer += int(m[1])
                    # 病毒亚种的kmer数计算包括该亚种及该亚种对应的种水平
                    if seqid2s1_tid.get(m[0]) == s1_tid or m[0] == s_tid:
                        spe_kmer += int(m[1])
            else:
                for mapping in ll[4].split(' '):
                    m = mapping.split(':')
                    total_kmer += int(m[1])
                    if seqid2speid.get(m[0]) == s_tid or m[0] == s_tid:
                        spe_kmer += int(m[1])
            rate = spe_kmer / total_kmer
            data.append(ll + [spe_kmer, rate])
            if kingdom == 'Viruses' and s1_tid:
                tax_id2rates[s1_tid].append(rate)
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
        df['kingdom'] = df.apply(lambda x: spe2kingdom.get(x[2]), axis=1)
        # get alignment for validating from viral
        df_viral = df[df['kingdom'] == 'Viruses']
        if not df_viral.empty:
# hy -start 20220516 23:34
            #for row in get_fail(df_viral, 's1_tid', int(top), seqid2s1_tid, spe_fail):
            #    wfp.write(row)
            if s1_tid: #hy
                for row in get_fail(df_viral, 's1_tid', int(top), seqid2s1_tid, spe_fail):
                    wfp.write(row)
            else:
                for row in get_fail(df_viral, 's_tid', int(top), seqid2speid, spe_fail): #hy
                    wfp.write(row)
# hy -end 20220516 23:34
        # get alignment for validating from non viral
        df_non_viral = df[df['kingdom'] != 'Viruses']
        if not df_non_viral.empty:
            for row in get_fail(df_non_viral, 's_tid', int(top), seqid2speid, spe_fail):
                wfp.write(row)


if __name__ == '__main__':
    main()
