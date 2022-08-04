#!/usr/bin/env python
# -* coding: utf-8 -*-

__author__ = 'sujiawei'
__email__ = 'su_jiawei163@163.com'
__date__ = '2022-04-14'
__version__ = 'V1.1'

import os
import re
import sys
import json
import click
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool, Process

import pandas as pd
import random
import warnings

Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg

warnings.filterwarnings('ignore')

KINGDOM2PREFIX = {
    'Bacteria': 'bacteria',
    'Eukaryota:Fungi': 'fungi',
    'SpecialPathogens': 'LY',
    'Eukaryota:Parasite': 'protozoa',
    'Eukaryota:Protozoa': 'protozoa',
    'Viruses': 'viral'
}


def get_taxid2info(db):
    # taxid info for viral
    taxid2info = {
        'seqtid2S1_tid': defaultdict(),
        'seqtid2S_tid': defaultdict(),
        'S1_tid2S_tid': defaultdict()
    }
    infile = Path(db).parent / 'seqTaxid2S1LevelTaxid'
    with open(infile, 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.strip().split('\t')
            if len(ll) >= 3:
                taxid2info['seqtid2S1_tid'][ll[0]] = ll[1]
                taxid2info['seqtid2S_tid'][ll[0]] = ll[2]
                taxid2info['S1_tid2S_tid'][ll[1]] = ll[2]
    return taxid2info


def species2taxid(db):
    spe2taxid = {}
    infile1 = Path(db).parent / 'all.Taxonomy.txt'
    infile2 = Path(db).parent / 'all.Taxonomy.other.txt'
    for file in infile1, infile2:
        df = pd.read_csv(file, sep='\t', dtype=str)
        for row in df.itertuples():
            if row[3] == 'Viruses':
                spe2taxid[row[10]] = row[1]
                spe2taxid[row[9]] = row[2]
            else:
                spe2taxid[row[9]] = row[2]

    infile3 = Path(db).parent / 'all.latin2chinese.INFO.txt'
    df = pd.read_csv(infile3, sep='\t', dtype=str)
    for row in df.itertuples():
        spe2taxid[row[3]] = row[2]

    return spe2taxid


def get_seqtid2stid(db):
    # taxid info for non viral
    infile = Path(db).parent / 'all.Taxonomy.txt'
    df = pd.read_csv(infile, sep='\t', dtype=str)
    seqtid2stid = {}
    for row in df[['Taxid', 'SpeciesTaxid']].itertuples():
        seqtid2stid[row[1]] = row[2]
        seqtid2stid[row[2]] = row[2]

    stid2group_tid = {}
    group_file = Path(db).parent / 'GroupList.info'
    df_grp = pd.read_csv(group_file, sep='\t', dtype=str)
    df_grp = df_grp[(df_grp['Focus'] == "*") | (df_grp['Focus'] == "+")]
    for row in df_grp[['GroupTaxid', 'SpeciesTaxid']].itertuples():
        stid2group_tid[row[1]] = row[1]
        stid2group_tid[row[2]] = row[1]

    return seqtid2stid, stid2group_tid


def get_seqtid2info(db):
    infile = Path(db) / f'{Path(db).name}.seqTaxid2Info.json'
    with open(infile, 'r') as rfp:
        seqtid2info = json.load(rfp)
    return seqtid2info


def custom_dict():
    ret_dict = defaultdict()
    for i in range(100):
        ret_dict[str(i + 1)] = 0
    return ret_dict


def get_cov(bam, tid2acc2cpn, kingdom, seqid2taxid, seqtid2info, seqtid2stid, stid2group_tid, df_bacteria):
    sps2cpn2depth = defaultdict(custom_dict)

    data = [line.split('\t') for line in os.popen(f'samtools view {bam}') if line.strip()]
    if not data:
        return kingdom, None
    df_bam = pd.DataFrame(data, columns=None)
    df_bam_sub = df_bam[[0, 2, 3]]
    df_bam_sub.columns = ['seq_id', 'ref_id', 'position']
    df_bam_sub['taxid'] = df_bam_sub.apply(lambda x: seqid2taxid.get(x['seq_id']), axis=1)
    df_bam_sub['S_tid'] = df_bam_sub.apply(lambda x: seqtid2stid.get(x['taxid']), axis=1)
    df_bam_sub['Grp_tid'] = df_bam_sub.apply(lambda x: stid2group_tid.get(x['taxid'], stid2group_tid.get(x['S_tid'])),
                                             axis=1)
    df_bam_sub = df_bam_sub.fillna('')

    # deal with species
    for stid in df_bam_sub['S_tid'].drop_duplicates():
        df_temp = df_bam_sub[df_bam_sub['S_tid'] == stid]
        for row in df_temp.itertuples():
            if stid in tid2acc2cpn:
                sps_name = seqtid2info[stid]['species']
                if row[2] in tid2acc2cpn[stid]:
                    for cpn, pos in tid2acc2cpn[stid][row[2]].items():
                        if pos[0] <= int(row[3]) <= pos[1]:
                            sps2cpn2depth[sps_name][cpn] += 1
                            break

    # deal with group/complex
    for grp_tid in df_bam_sub['Grp_tid'].drop_duplicates():
        grp_rc = len(df_bacteria[df_bacteria['tax_id'] == grp_tid])
        df_temp = df_bam_sub[df_bam_sub['Grp_tid'] == grp_tid]
        df_temp = df_temp.iloc[:grp_rc, :]
        for row in df_temp.itertuples():
            if grp_tid in tid2acc2cpn:
                sps_name = seqtid2info[grp_tid]['species']
                if row[2] in tid2acc2cpn[grp_tid]:
                    for cpn, pos in tid2acc2cpn[grp_tid][row[2]].items():
                        if pos[0] <= int(row[3]) <= pos[1]:
                            sps2cpn2depth[sps_name][cpn] += 1
                            break

    return kingdom, sps2cpn2depth


def get_cov4sub_virus(bam, tid2acc2cpn, kingdom, seqid2taxid, df_virus, taxid2info):
    """ 根据病毒的检出结果及分型相关的注释表，统计病毒的覆盖度 """
    data = [line.split('\t') for line in os.popen(f'samtools view {bam}') if line.strip()]
    if not data:
        return kingdom, None
    df_bam = pd.DataFrame(data, columns=None)
    df_bam_sub = df_bam[[0, 2, 3]]
    df_bam_sub.columns = ['seq_id', 'ref_id', 'position']
    df_bam_sub['taxid'] = df_bam_sub.apply(lambda x: seqid2taxid.get(x['seq_id']), axis=1)
    df_bam_sub['S_tid'] = df_bam_sub.apply(lambda x: taxid2info['seqtid2S_tid'].get(x['taxid']), axis=1)
    df_bam_sub['S1_tid'] = df_bam_sub.apply(lambda x: taxid2info['seqtid2S1_tid'].get(x['taxid']), axis=1)
    df_bam_sub['refid2S_tid'] = df_bam_sub.apply(lambda x: taxid2info['seqtid2S1_tid'].get(x['ref_id'].split('|')[0]),
                                                 axis=1)
    sps2cpn2depth = defaultdict(custom_dict)

    df_virus = df_virus[['tax_id', 'S_tid', 'SpeciesReads', 'SpeciesSN']]
    for row in df_virus.itertuples():
        # species level
        if row[1] == row[2]:
            df_temp = df_bam_sub[(df_bam_sub['S_tid'] == row[2]) | (df_bam_sub['taxid'] == row[2])]
            for bam_row in df_temp.itertuples():
                if row[2] in tid2acc2cpn:
                    if bam_row[2] in tid2acc2cpn[row[2]]:
                        for cpn, pos in tid2acc2cpn[row[2]][bam_row[2]].items():
                            # 若比对位置位于基因组某一组分的区间内，则该组分的覆盖深度+1
                            if pos[0] <= int(bam_row[3]) <= pos[1]:
                                sps2cpn2depth[row[4]][cpn] += 1
                                break
        # subspecies level
        else:
            df_temp_subsps = df_bam_sub[df_bam_sub['S1_tid'] == row[1]]
            s1_reads_cnt = len(df_temp_subsps)
            if s1_reads_cnt < int(row[3]):
                diff_val = int(row[3]) - s1_reads_cnt
                df_temp_sps = df_bam_sub[(df_bam_sub['taxid'] == row[2])]
                idx_temp = list(df_temp_sps.index)
                random.shuffle(idx_temp)
                for idx in idx_temp[:diff_val]:
                    # 将差值部分的比对结果由种改为亚种
                    df_bam_sub.loc[idx, 'taxid'] = row[1]
                    df_bam_sub.loc[idx, 'S1_tid'] = row[1]
                df_temp_subsps = df_bam_sub[df_bam_sub['S1_tid'] == row[1]]
            else:
                df_temp_subsps = df_temp_subsps.iloc[:s1_reads_cnt]
            for bam_row in df_temp_subsps.itertuples():
                if bam_row[7] in tid2acc2cpn:
                    tid_idx = bam_row[7]
                elif row[1] in tid2acc2cpn:
                    tid_idx = row[1]
                elif row[2] in tid2acc2cpn:
                    tid_idx = row[2]
                else:
                    continue
                if bam_row[2] in tid2acc2cpn[tid_idx]:
                    for cpn, pos in tid2acc2cpn[tid_idx][bam_row[2]].items():
                        if pos[0] <= int(bam_row[3]) <= pos[1]:
                            sps2cpn2depth[row[4]][cpn] += 1
                            break
    return kingdom, sps2cpn2depth


def stat_cov(indir, db, spe2taxid, seqtid2info, seqtid2stid, stid2group_tid, parallel, name, outdir):
    # 获取病毒亚种与种水平taixd关系表
    taxid2info = get_taxid2info(db)
    # 获取物种检出结果总表
    df_detect = pd.read_excel(Path(indir).parent / 'Total_Detail.xlsx', dtype=str)
    # 获取病毒检出结果表
    df_virus = df_detect[df_detect['Kingdom'] == 'Viruses']
    if not df_virus.empty:
        df_virus['tax_id'] = df_virus.apply(lambda x: spe2taxid.get(x['SpeciesSN']), axis=1)
        df_virus['S_tid'] = df_virus.apply(lambda x: taxid2info['S1_tid2S_tid'].get(x['tax_id'], x['tax_id']), axis=1)

    # 获取细菌检出结果表
    df_bacteria = df_detect[df_detect['Kingdom'] == 'Bacteria']
    df_bacteria['tax_id'] = df_bacteria.apply(lambda x: spe2taxid.get(x['SpeciesSN']), axis=1)

    # 获取序列比对到的参考基因组的taxid
    seqid2taxid = {}
    with open(Path(indir).parent / f'{name}.unique.kout', 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.split('\t')
            seqid2taxid[ll[1]] = ll[2]

    # 进程池，规定运行时并行子进程的数量
    pool = Pool(processes=int(parallel))
    # 存放函数的返回值，即子进程的返回值
    results = []
    for bam in Path(indir).glob('*/*bam'):
        gtid = bam.name.split('.')[0]
        kingdom = bam.parent.name
        if gtid == 'NULL':
            gtid = 'NULL_' + kingdom
        if bam.stat().st_size == 0:
            continue
        # 获取genus_taxid，根据genus_taxid获取注释信息
        tid2acc2cpn_jsn = Path(db) / gtid / f'{gtid}.tid2acc2cpn.json'
        tid2acc2cpn_fp = tid2acc2cpn_jsn.open('r')
        tid2acc2cpn = json.load(tid2acc2cpn_fp)
        tid2acc2cpn_fp.close()
        # 并行处理任务
        if kingdom == 'Viruses':
            results.append(pool.apply_async(get_cov4sub_virus,
                                            args=(bam, tid2acc2cpn, kingdom, seqid2taxid, df_virus, taxid2info)))
        else:
            results.append(pool.apply_async(get_cov,
                                            args=(bam, tid2acc2cpn, kingdom, seqid2taxid, seqtid2info, seqtid2stid,
                                                  stid2group_tid, df_bacteria)))
    pool.close()
    pool.join()

    kingdom2wfp = defaultdict()
    for k, v in KINGDOM2PREFIX.items():
        kingdom2wfp[k] = Path(outdir).joinpath(f'{v}.{name}.cov').open('w')

    # 将子进程返回值按指定格式输出
    for ret in results:
        kingdom, sps2cpn2depth = ret.get()
        if not sps2cpn2depth:
            continue
        for sps_name, cpn2depth in sps2cpn2depth.items():
            for cpn, depth in cpn2depth.items():
                kingdom2wfp[kingdom].write(f'{sps_name}\t{str(cpn)}\t{str(depth)}\n')

    for wfp in kingdom2wfp.values():
        wfp.close()

    for cov in Path(outdir).glob('*cov'):
        if cov.stat().st_size == 0:
            os.system(f'rm {str(cov)}')


def draw_cov_plot(dir, script):
    for cov in Path(dir).glob('*cov'):
        kingdom = cov.name.split('.')[0]
        outdir = Path(dir) / f'{kingdom}_readscount_coverage_plot'
        if not outdir.exists():
            os.system(f'mkdir -p {str(outdir)}')

        os.system(f'{cfg.Rscript} {script} {str(cov)} {str(outdir)}')

        df_cov = pd.read_csv(cov, sep='\t', dtype=str, header=None)
        for sps_name in set(df_cov[0].tolist()):
            sps_name = sps_name.replace(' ', '_')
            if '(' in sps_name or ')' in sps_name:
                sps_name = sps_name.replace('/', '_')
                src_file = outdir / (sps_name.replace('(', '\\(').replace(')', '\\)') + '.coverage.png')
                tgt_file = outdir / (sps_name + '.coverage.png')
                os.system(f'cp \'{str(src_file)}\' \'{str(tgt_file)}\'')


@click.command(no_args_is_help=True)
@click.option('-i', '--indir', help='Bwa比对结果文件所在目录, eg: BwaAlign/')
@click.option('--db', help='absolute path to bwa database')
@click.option('-p', '--parallel', help='num of parallel jobs', default=10)
@click.option('-s', '--script', help='Script for drawing coverage plots')
@click.option('-o', '--outdir', help='output dir, eg: CovPlots/')
def main(indir, db, parallel, script, outdir):
    """ 根据指定路径下的bam文件，统计各物种对应的参考基因组（共分为100组分）上的覆盖深度 """
    # 获取BWA数据库中seq_taxid对应的species_taxid, genus_taxid, species_name
    seqtid2info = get_seqtid2info(db)

    # 获取PIDB数据库中的物种分类信息
    spe2taxid = species2taxid(db)
    seqtid2stid, stid2group_tid = get_seqtid2stid(db)

    if not Path(outdir).exists():
        os.system(f'mkdir -p {str(outdir)}')

    # 统计覆盖度
    name = Path(indir).parent.parent.name
    stat_cov(indir, db, spe2taxid, seqtid2info, seqtid2stid, stid2group_tid, parallel, name, outdir)

    # 画图
    draw_cov_plot(outdir, script)


if __name__ == '__main__':
    main()
