#!/usr/bin/env python

# ===============================================================================
#   FileName: QC.stat.pl
#   Author  : qiuwancen
#   E-Mail  : 972538446@qq.com
#   Version : 1.0
#   Date    : 2022-02-28
# ===============================================================================

import os
import re
import sys
import json
import gzip
import click
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict


def check_file(path):
    if not Path(path).exists():
        raise FileExistsError(f"{str(path)} not exists!")


def custom_dict():
    return {
        'Sample': '',
        'Type': '',
        'Raw Reads (M)': '',
        'Raw Bases': '',
        'Raw Bases Q20 (%)': '',
        'Raw Bases Q30 (%)': '',
        'Clean Reads': '',
        'Clean Bases': '',
        'Clean Bases Q20 (%)': '',
        'Clean Bases Q30 (%)': '',
        'Duplication Rate (%)': '',
        'Reads with Adapter': '',
        'Bases with Adapter': '',
        'Reads with N': '',
        'Reads with Low Quality': '',
        'Read with Low Complexity': '',
        'Reads filtered due to too short': '',
        'Host Reads': '',
        'Host Rate (%)': '',
        'Internal Control Reads': '',
        'Internal Control Rate (%)': '',
        'Internal Control Gene18 Average Depth': '',
        'Internal Control Gene18 Coverage (%)': '',
        'Internal Control Gene43 Average Depth': '',
        'Internal Control Gene43 Coverage (%)': '',
        'Classified_Reads_Number': '',
        'Classified_Rate (%)': '',
        'Unclassified_Reads_Number': '',
        'Unclassified_Rate (%)': '',
        'Viruses_Reads_Number': '',
        'Viruses_Rate (%)': '',
        'Bacteria_Reads_Number': '',
        'Bacteria_Rate (%)': '',
        'Fungi_Reads_Number': '',
        'Fungi_Rate (%)': '',
        'Protozoa_Reads_Number': '',
        'Protozoa_Rate (%)': '',
        'Metazoa_Parasite_Reads_Number': '',
        'Metazoa_Parasite_Rate (%)': '',
    }


def reformat(obj):
    if isinstance(obj, float):
        return round(obj, 6)
    else:
        return obj


@click.command(no_args_is_help=True)
@click.option('-w', '--workdir', help='分析结果路径')
@click.option('-s', '--samples', help='样本列表')
def main(workdir, samples):
    qc_dir = Path(workdir) / '01.QC/01.Fastp'
    rh_dir = Path(workdir) / '01.QC/02.RemoveHost'
    tx_dir = Path(workdir) / '02.Annotation/01.Taxonomy'
    check_file(samples)
    with open(samples, 'r') as rfp:
        sample2library = {}
        for line in rfp.readlines():
            if line.strip():
                ll = line.split('\t')
                sample2library[ll[0]] = ll[1]

    qc_stat = defaultdict(custom_dict)
    qc2list = defaultdict(list)

    for sample, library in sample2library.items():
        qc_stat[sample]['Sample'] = sample
        qc_stat[sample]['Type'] = library
        # fastp QC
        jsn_file = qc_dir / f'{sample}/{library}/{sample}.json'
        check_file(jsn_file)
        fastp_info = json.load(jsn_file.open('r'))
        qc_stat[sample]['Raw Reads (M)'] = fastp_info["summary"]["before_filtering"]["total_reads"] / 10 ** 6
        qc_stat[sample]['Raw Bases'] = fastp_info["summary"]["before_filtering"]["total_bases"]
        qc_stat[sample]['Raw Bases Q20 (%)'] = fastp_info["summary"]["before_filtering"]["q20_rate"] * 100
        qc_stat[sample]['Raw Bases Q30 (%)'] = fastp_info["summary"]["before_filtering"]["q30_rate"] * 100
        qc_stat[sample]['Clean Reads'] = fastp_info["summary"]["after_filtering"]["total_reads"]
        qc_stat[sample]['Clean Bases'] = fastp_info["summary"]["after_filtering"]["total_bases"]
        qc_stat[sample]['Clean Bases Q20 (%)'] = fastp_info["summary"]["after_filtering"]["q20_rate"] * 100
        qc_stat[sample]['Clean Bases Q30 (%)'] = fastp_info["summary"]["after_filtering"]["q30_rate"] * 100
        qc_stat[sample]['Duplication Rate (%)'] = fastp_info["duplication"]["rate"] * 100
        qc_stat[sample]['Reads with Adapter'] = fastp_info["adapter_cutting"]["adapter_trimmed_reads"]
        qc_stat[sample]['Bases with Adapter'] = fastp_info["adapter_cutting"]["adapter_trimmed_bases"]
        qc_stat[sample]['Reads with N'] = fastp_info["filtering_result"]["too_many_N_reads"]
        qc_stat[sample]['Reads with Low Quality'] = fastp_info["filtering_result"]["low_quality_reads"]
        qc_stat[sample]['Read with Low Complexity'] = fastp_info["filtering_result"]["low_complexity_reads"]
        qc_stat[sample]['Reads filtered due to too short'] = fastp_info["filtering_result"]["too_short_reads"]
        # batch stat
        qc2list['RawReads'].append(fastp_info["summary"]["before_filtering"]["total_reads"])
        qc2list['CleanReads'].append(fastp_info["summary"]["after_filtering"]["total_reads"])
        qc2list['CleanQ30'].append(fastp_info["summary"]["after_filtering"]["q30_rate"])
        # Remove host stat
        rh_file = rh_dir / sample / library / f'{sample}.Homo.log'
        check_file(rh_file)
        with open(rh_file, 'r') as rfp:
            for line in rfp.readlines():
                if re.match('\d+ \+ \d+ mapped ', line):
                    ll = line.split(' ')
                    qc_stat[sample]['Host Reads'] = int(ll[0])
                    break
        # Internal control
        ic_file = rh_dir / sample / library / f'{sample}.delIT.log'
        if ic_file.exists():
            with open(ic_file, 'r') as rfp:
                for line in rfp.readlines():
                    if re.match('\d+ \+ \d+ mapped ', line):
                        ll = line.split(' ')
                        qc_stat[sample]['Internal Control Reads'] = ll[0]
                        qc_stat[sample]['Internal Control Rate (%)'] = ll[4].replace('(', '', 1)
                        break
        # Internal control coverage and depth
        rg_file = rh_dir / sample / library / 'region.tsv.gz'
        if rg_file.exists():
            with gzip.open(rg_file, 'rb') as rfp:
                for line in rfp:
                    ll = str(line, 'utf-8').split('\t')
                    if ll[1] == '2297':
                        qc_stat[sample]['Internal Control Gene18 Average Depth'] = ll[3]
                        qc_stat[sample]['Internal Control Gene18 Coverage (%)'] = ll[5]
                    if ll[1] == '4238':
                        qc_stat[sample]['Internal Control Gene43 Average Depth'] = ll[3]
                        qc_stat[sample]['Internal Control Gene43 Coverage (%)'] = ll[5]
        # Taxonomy stat
        ts_file = tx_dir / sample / library / 'Taxonomy_Summary.txt'
        check_file(ts_file)
        df = pd.read_csv(ts_file, sep='\t', dtype=str)
        clean_reads = qc_stat[sample]['Clean Reads']
        qc_stat[sample]['Classified_Reads_Number'] = int(df.loc[0, 'Classified_Reads_Number']) - int(df.loc[0, 'Human'])
        qc_stat[sample]['Classified_Rate (%)'] = qc_stat[sample]['Classified_Reads_Number'] / clean_reads * 100
        qc_stat[sample]['Unclassified_Reads_Number'] = int(df.loc[0, 'Unclassified_Reads_Number'])
        qc_stat[sample]['Unclassified_Rate (%)'] = qc_stat[sample]['Unclassified_Reads_Number'] / clean_reads * 100
        qc_stat[sample]['Host Reads'] += int(df.loc[0, 'Human'])
        qc_stat[sample]['Host Rate (%)'] = qc_stat[sample]['Host Reads'] / clean_reads * 100
        qc_stat[sample]['Viruses_Reads_Number'] = int(df.loc[0, 'Viruses'])
        qc_stat[sample]['Viruses_Rate (%)'] = qc_stat[sample]['Viruses_Reads_Number'] / clean_reads * 100
        qc_stat[sample]['Bacteria_Reads_Number'] = int(df.loc[0, 'Bacteria'])
        qc_stat[sample]['Bacteria_Rate (%)'] = qc_stat[sample]['Bacteria_Reads_Number'] / clean_reads * 100
        qc_stat[sample]['Fungi_Reads_Number'] = int(df.loc[0, 'Fungi'])
        qc_stat[sample]['Fungi_Rate (%)'] = qc_stat[sample]['Fungi_Reads_Number'] / clean_reads * 100
        qc_stat[sample]['Protozoa_Reads_Number'] = int(df.loc[0, 'Protozoa'])
        qc_stat[sample]['Protozoa_Rate (%)'] = qc_stat[sample]['Protozoa_Reads_Number'] / clean_reads * 100
        qc_stat[sample]['Metazoa_Parasite_Reads_Number'] = int(df.loc[0, 'Metazoa_Parasite'])
        qc_stat[sample]['Metazoa_Parasite_Rate (%)'] = qc_stat[sample]['Metazoa_Parasite_Reads_Number'] / clean_reads * 100

    df_out = pd.DataFrame().from_dict(qc_stat)
    for col in df_out.columns:
        df_out[col] = df_out.apply(lambda x: reformat(x[col]), axis=1)

    # output qc for each sample
    for sample in sample2library.keys():
        sample_dir = Path(workdir) / 'plugin_out/Result' / sample
        if not sample_dir.exists():
            os.system(f'mkdir -p {str(sample_dir)}')
        sample_out_tsv = sample_dir / 'QC.Stat'
        df_out[[sample]].to_csv(sample_out_tsv, sep='\t', header=False)
        sample_out_xlsx = sample_dir / 'QC.Stat.xlsx'
        df_out[[sample]].to_excel(sample_out_xlsx, header=False)

    # stat sequencing info
    
    batch_qc = {'RunID': 'Test', 'Number of samples': len(sample2library)}
    df_smp = pd.read_csv(samples, sep='\t', dtype=str)
    smr_tbl = Path(df_smp.iloc[0, 2]).parent / 'summaryTable.csv'
    df_smr = pd.read_csv(smr_tbl, dtype=str)
    batch_qc['ChipProductivity (%)'] = float(df_smr.iloc[4, 1])
    batch_qc['SplitRate (%)'] = float(df_smr.iloc[8, 1])
    batch_qc['Raw Q30 (%)'] = float(df_smr.iloc[7, 1])
    batch_qc['TotalRawReads'] = sum(qc2list['RawReads'])
    batch_qc['MaxRawReads'] = max(qc2list['RawReads'])
    batch_qc['MinRawReads'] = min(qc2list['RawReads'])
    batch_qc['MeanRawReads'] = np.mean(qc2list['RawReads'])
    batch_qc['TotalCleanReads'] = sum(qc2list['CleanReads'])
    batch_qc['MaxCleanReads'] = max(qc2list['CleanReads'])
    batch_qc['MinCleanReads'] = min(qc2list['CleanReads'])
    batch_qc['MeanCleanReads'] = np.mean(qc2list['CleanReads'])
    batch_qc['Clean Q30 (%)'] = np.mean(qc2list['CleanQ30']) * 100
    batch_qc['Effective (%)'] = batch_qc['MaxCleanReads'] / batch_qc['TotalRawReads'] * 100
    batch_qc['Number of QC-PASS samples'] = len(sample2library)
    batch_qc['Number of QC-FAIL samples'] = len(sample2library)

    # output
    out_dir = Path(workdir) / 'plugin_out/Result'
    if not out_dir.exists():
        os.system(f'mkdir -p {str(out_dir)}')
    out_tsv = out_dir / 'QC.Stat'
    df_out.to_csv(out_tsv, sep='\t', header=False)

    # merge df
    with pd.ExcelWriter(out_dir / 'QC.Stat.xlsx') as writer:
        df_batch = pd.DataFrame.from_dict(data=batch_qc, orient='index')
        df_batch.to_excel(writer, sheet_name='Sheet1', header=False)
        df_out.to_excel(writer, sheet_name='Sheet2', header=False)


if __name__ == '__main__':
    main()
