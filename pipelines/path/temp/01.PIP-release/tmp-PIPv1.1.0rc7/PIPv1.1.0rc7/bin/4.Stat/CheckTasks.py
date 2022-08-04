#!/usr/bin/python
import os
import sys
import json
import time
import click
import pandas as pd
from pathlib import Path
from collections import defaultdict


############################ developer info ##############################
__AUTHOR__ = "Sujiawei"
__EMAIL__ = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__ = "2022-06-08"
##########################################################################


def get_sample(smp2lib, sh_file):
    """ 根据task.sh中的命令行定位样本编号 """
    with open(sh_file) as rfp:
        cmd_list = rfp.readlines()
    for smp in smp2lib.keys():
        for line in cmd_list:
            if smp in line:
                return smp


def get_task2status(infile):
    """ 获取任务脚本的分析状态 """
    task2status = defaultdict(set)
    status_set = {'Start', 'Error', 'Finish'}
    with open(infile, 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.strip().split(': ')
            if ll[1] in status_set:
                if ll[1] == 'Error':
                    task2status[ll[3]].add(ll[1])
                else:
                    task2status[ll[2]].add(ll[1])
    return task2status


def get_error(infile):
    """ 获取标准错误输出的日志数据 """
    infile = Path(infile)
    for err_file in infile.parent.glob(infile.name + '*.e'):
        if err_file:
            with open(err_file, 'r') as rfp:
                return ''.join(rfp.readlines())
        else:
            return 'Error file not found.'


def get_note(task2status, smp2lib):
    """ 获取样本对应的分析状态 """
    label = 'Info'
    log_text = ''
    for task, status in task2status.items():
        if 'Error' in status:
            note = '分析异常'
            label = 'Error'
            log_text = get_error(task)
        elif 'Finish' in status:
            note = '分析完成'
        else:
            note = '分析中'
        smp = get_sample(smp2lib, task)
        yield task, smp, note, label, log_text


def check_fastp(fastp_dir, sample, library):
    """ 检查质控结果文件是否齐全并返回信息 """
    gz_file = Path(fastp_dir) / sample / library / f'{sample}.fq.gz'
    html_file = Path(fastp_dir) / sample / library / f'{sample}.html'
    json_file = Path(fastp_dir) / sample / library / f'{sample}.json'
    log_file = Path(fastp_dir) / sample / library / f'{sample}.log'
    if not gz_file.exists() or not html_file.exists() or not json_file.exists():
        if log_file.exists():
            with open(log_file, 'r') as rfp:
                return ''.join(rfp.readlines())
        else:
            return f'{str(log_file)} not exists.'
    else:
        return 'No error...'


def check_tax(tax_dir, sample, library):
    """ 检查物种鉴定结果文件是否齐全并返回信息 """
    file_path = Path(tax_dir) / sample / library / 'Total_Detail.txt'
    if not file_path.exists():
        return 'Error', f'检出结果表不存在: {str(file_path)}'
    else:
        df = pd.read_csv(file_path, sep='\t', dtype=str)
        if not df.empty:
            gz_lst = [gz for gz in file_path.parent.glob('*result.fq.gz')]
            if gz_lst:
                return 'Info', ''
            else:
                return 'Error', f'检出结果序列不存在: {str(file_path.parent)}/*result.fq.gz'


def check_cov(tax_dir, sample, library):
    """ 检查物种覆盖图结果文件是否齐全并返回信息 """
    plot_dir = Path(tax_dir) / sample / library / 'Coverage/04.plot'
    is_err = True
    for png in plot_dir.glob('*/*png'):
        if png.exists():
            is_err = False
            return is_err
    return is_err


def check_result(sample, smp_dir, workdir):
    """ 检查统计结果文件是否齐全并返回信息 """
    log = []
    file_lst = ['QC.Stat.xlsx', 'raw_summary.xls', 'Taxonomy_Summary.txt',
                'ARG.cov.pick.xlsx', 'ARG.cov.pick.txt', 'VF.xlsx', 'VF.txt',
                'Pathogeny_Detail.txt', 'Pathogeny_Detail.xlsx',
                'Total_Detail.txt', 'Total_Detail.xlsx', f'{sample}.merge.xls']
    for file in file_lst:
        file_path = Path(smp_dir) / file
        if not file_path.exists():
            log.append(f'Error: {str(file_path)} 不存在! 请检查上游分析结果或统计步骤: Stat.py')

    plot_dirs = [d for d in os.listdir(smp_dir) if 'plot' in d]
    if len(plot_dirs) == 0:
        log.append('Error: 无覆盖图文件夹, Coverage.py或Stat.py步骤可能存在异常，请检查!')
    else:
        has_png = False
        for plot_dir in plot_dirs:
            if len(os.listdir(smp_dir / plot_dir)) > 0:
                has_png = True
                break
        if not has_png:
            log.append(f'Error: 无覆盖图文件, Coverage.py或Stat.py步骤可能存在异常, 请检查!')

    stat_log = Path(workdir) / '04.Stat.log'
    if stat_log.exists():
        with open(stat_log, 'r') as rfp:
            for line in rfp.readlines():
                if line.startswith('Warning: '):
                    log.append(line)
    else:
        log.append(f'Error: {str(stat_log)} 统计步骤日志文件不存在, 请检查!')

    log_text = '\n'.join(log)
    if 'Error: ' in log_text:
        return 'Error', log_text
    elif 'Warning: ' in log_text:
        return 'Warning', log_text
    else:
        return 'Info', log_text


def custom_dict():
    return {
        '00.Merge': {'comment': '数据合并', 'note': '-', 'label': 'Info', 'log': ''},
        '01.Fastp': {'comment': 'QC质控', 'note': '-', 'label': 'Info', 'log': ''},
        '02.RemoveHost': {'comment': '去宿主', 'note': '-', 'label': 'Info', 'log': ''},
        '03.Taxonomy-Classify': {'comment': '物种鉴定-比对统计', 'note': '-', 'label': 'Info', 'log': ''},
        '04.Taxonomy-Coverage': {'comment': '物种鉴定-覆盖度', 'note': '-', 'label': 'Info', 'log': ''},
        '05.ResistantGene': {'comment': '耐药检测', 'note': '-', 'label': 'Info', 'log': ''},
        '06.VF': {'comment': '毒力分析', 'note': '-', 'label': 'Info', 'log': ''},
        '07.Stat': {'comment': '批次统计', 'note': '-', 'label': 'Info', 'log': ''}
    }


def workflow(samples, workdir, tag):
    """ 检查各个步骤的分析状态，反馈分析进度 """
    workdir = Path(workdir)
    merge_dir = workdir / '00.Merge'
    fastp_dir = workdir / '01.QC/01.Fastp'
    rm_host_dir = workdir / '01.QC/02.RemoveHost'  # 去除宿主
    tax_dir = workdir / '02.Annotation/01.Taxonomy'  # 物种鉴定
    arg_dir = workdir / '02.Annotation/02.ResistanceGene'  # 耐药基因
    vf_dir = workdir / '02.Annotation/03.VFDB'  # 毒力分析
    res_dir = workdir / 'plugin_out/Result'  # 结果目录

    df_smp = pd.read_csv(samples, sep='\t', dtype=str, header=None)
    smp2lib = {row[1]: row[2] for row in df_smp.itertuples()}

    result = {
        'batchSchedule': {'comment': '批次分析进度', 'note': '-', 'log': '',
                          'finished': 0, 'failed': 0, 'processing': 0, 'total': len(smp2lib),
                          'finished_samples': [], 'failed_samples': [], 'processing_samples': [],
                          'batch_id': Path(workdir).name, 'batch_path': str(workdir)},
        'sampleSchedule': defaultdict(custom_dict),
        'label': 'Info'
    }

    # deal with 00.Merge
    for mg_log in (merge_dir / 'shell').glob('Merge.sh.*.log'):
        if mg_log.exists():
            task2status = get_task2status(mg_log)
            for task, smp, note, label, log_text in get_note(task2status, smp2lib):
                result['sampleSchedule'][smp]['00.Merge']['note'] = note
                result['sampleSchedule'][smp]['00.Merge']['label'] = label
                result['sampleSchedule'][smp]['00.Merge']['log'] = log_text

    # deal with 01.QC/01.Fastp
    for fsp_log in (fastp_dir / 'shell').glob('Fastp.sh.*.log'):
        if fsp_log.exists():
            task2status = get_task2status(fsp_log)
            for task, status in task2status.items():
                smp = get_sample(smp2lib, task)
                label = 'Info'
                if 'Error' in status:
                    log_text = get_error(task)
                    if not log_text:
                        log_text = check_fastp(fastp_dir, smp, smp2lib[smp])
                    note = '分析异常'
                    label = 'Error'
                elif 'Finish' in status:
                    note = '分析完成'
                else:
                    note = '分析中'
                result['sampleSchedule'][smp]['01.Fastp']['note'] = note
                result['sampleSchedule'][smp]['01.Fastp']['label'] = label
                result['sampleSchedule'][smp]['01.Fastp']['log'] = log_text

    # deal with 01.QC/02.RemoveHost
    for rh_log in (rm_host_dir / 'shell').glob('RemoveHost.sh.*.log'):
        if rh_log.exists():
            task2status = get_task2status(rh_log)
            for task, smp, note, label, log_text in get_note(task2status, smp2lib):
                result['sampleSchedule'][smp]['02.RemoveHost']['note'] = note
                result['sampleSchedule'][smp]['02.RemoveHost']['label'] = label
                result['sampleSchedule'][smp]['02.RemoveHost']['log'] = log_text
                if note == '分析完成':
                    rh_gz = rm_host_dir / smp / smp2lib[smp] / f'{smp}.fq.gz'
                    if not rh_gz.exists():
                        result['sampleSchedule'][smp]['02.RemoveHost']['note'] = '分析异常'
                        result['sampleSchedule'][smp]['02.RemoveHost']['label'] = 'Error'
                        err_log_file = rm_host_dir / smp / smp2lib[smp] / f'{smp}.SNAP.log'
                        with open(err_log_file, 'r') as rfp:
                            err_log = ''.join(rfp.readlines())
                        result['sampleSchedule'][smp]['02.RemoveHost']['log'] = get_error(err_log)

    # deal with 02.Annotation/01.Taxonomy
    for krk_log in (tax_dir / 'shell').glob('Kraken.sh.*.log'):
        if krk_log.exists():
            task2status = get_task2status(krk_log)
            for task, smp, note, label, log_text in get_note(task2status, smp2lib):
                result['sampleSchedule'][smp]['03.Taxonomy-Classify']['note'] = note
                result['sampleSchedule'][smp]['03.Taxonomy-Classify']['label'] = label
                result['sampleSchedule'][smp]['03.Taxonomy-Classify']['log'] = log_text
                if note == '分析完成':
                    label, log_text = check_tax(tax_dir, smp, smp2lib[smp])
                    if label == 'Error':
                        result['sampleSchedule'][smp]['03.Taxonomy-Classify']['note'] = '分析异常'
                        result['sampleSchedule'][smp]['03.Taxonomy-Classify']['label'] = label
                        err_log = '\n\n'.join([log_text, get_error(Path(task))])
                        result['sampleSchedule'][smp]['03.Taxonomy-Classify']['log'] = err_log
    for cov_log in (tax_dir / 'shell').glob('Coverage.sh.*.log'):
        if cov_log.exists():
            task2status = get_task2status(cov_log)
            for task, smp, note, label, log_text in get_note(task2status, smp2lib):
                result['sampleSchedule'][smp]['04.Taxonomy-Coverage']['note'] = note
                result['sampleSchedule'][smp]['04.Taxonomy-Coverage']['label'] = label
                result['sampleSchedule'][smp]['04.Taxonomy-Coverage']['log'] = log_text
                if note == '分析完成':
                    is_err = check_cov(tax_dir, smp, smp2lib[smp])
                    if is_err:
                        result['sampleSchedule'][smp]['04.Taxonomy-Coverage']['note'] = '分析异常'
                        result['sampleSchedule'][smp]['04.Taxonomy-Coverage']['label'] = 'Error'
                        result['sampleSchedule'][smp]['04.Taxonomy-Coverage']['log'] = get_error(Path(task))

    # deal with 02.Annotation/02.ResistanceGene
    for arg_log in (arg_dir / 'shell').glob('ResistanceGene.sh.*.log'):
        if arg_log.exists():
            task2status = get_task2status(arg_log)
            for task, smp, note, label, log_text in get_note(task2status, smp2lib):
                result['sampleSchedule'][smp]['05.ResistantGene']['note'] = note
                result['sampleSchedule'][smp]['05.ResistantGene']['label'] = label
                result['sampleSchedule'][smp]['05.ResistantGene']['log'] = log_text
                if note == '分析完成':
                    arg_res_file = arg_dir / smp / smp2lib[smp] / 'ARG.cov.pick.txt'
                    if not arg_res_file.exists():
                        result['sampleSchedule'][smp]['05.ResistantGene']['note'] = '分析异常'
                        result['sampleSchedule'][smp]['05.ResistantGene']['label'] = 'Error'
                        result['sampleSchedule'][smp]['05.ResistantGene']['log'] = get_error(Path(task))

    # deal with 02.Annotation/03.VFDB
    for vf_log in (vf_dir / 'shell').glob('VFDB.sh.*.log'):
        if vf_log.exists():
            task2status = get_task2status(vf_log)
            for task, smp, note, label, log_text in get_note(task2status, smp2lib):
                result['sampleSchedule'][smp]['06.VF']['note'] = note
                result['sampleSchedule'][smp]['06.VF']['label'] = label
                result['sampleSchedule'][smp]['06.VF']['log'] = log_text
                if note == '分析完成':
                    vf_res_file = vf_dir / smp / smp2lib[smp] / 'VF.txt'
                    if not vf_res_file.exists():
                        result['sampleSchedule'][smp]['06.VF']['note'] = '分析异常'
                        result['sampleSchedule'][smp]['06.VF']['label'] = 'Error'
                        result['sampleSchedule'][smp]['06.VF']['log'] = get_error(Path(task))

    # deal with plugin_out/Result
    successful_smp = set()
    processing_smp = set()
    failed_smp = set()
    for analysis in custom_dict().keys():
        if analysis == '07.Stat':
            continue
        ana_failed_smps = set()
        ana_success_smps = set()
        for smp in result['sampleSchedule'].keys():
            if result['sampleSchedule'][smp][analysis]['note'] == '分析异常':
                failed_smp.add(smp)
                ana_failed_smps.add(smp)
            elif result['sampleSchedule'][smp][analysis]['note'] == '分析中':
                processing_smp.add(smp)
            else:
                ana_success_smps.add(smp)
        if len(ana_failed_smps):
            result['label'] = 'Error'
            result['batchSchedule']['note'] = '分析异常'
            old = result['batchSchedule']['log']
            new = f'步骤: {analysis} 分析异常, 请排查此步骤!'
            result['batchSchedule']['log'] = '\n'.join([old, new])
        elif len(ana_success_smps) == len(smp2lib):
            old = result['batchSchedule']['log']
            new = f'步骤: {analysis} 分析完成, 所有样本均分析成功.'
            result['batchSchedule']['log'] = '\n'.join([old, new])
        else:
            result['batchSchedule']['processing'] = len(processing_smp)
            old = result['batchSchedule']['log']
            new = f'步骤: {analysis} 分析中...'
            result['batchSchedule']['log'] = '\n'.join([old, new])

    if tag:  # by2.done is existed
        # sample stat
        for smp in result['sampleSchedule'].keys():
            smp_dir = res_dir / smp
            if not smp_dir.exists():
                result['sampleSchedule'][smp]['07.Stat']['note'] = '分析异常'
            else:
                label, log_text = check_result(smp, smp_dir, workdir)
                if label != 'Error':
                    successful_smp.add(smp)
                    result['sampleSchedule'][smp]['07.Stat']['note'] = '分析完成'
                else:
                    failed_smp.add(smp)
                    result['sampleSchedule'][smp]['07.Stat']['note'] = '分析异常'
                result['sampleSchedule'][smp]['07.Stat']['label'] = label
                result['sampleSchedule'][smp]['07.Stat']['log'] = log_text
        # batch stat
        if not res_dir.exists():
            result["batchSchedule"]['note'] = '分析异常'
        if len(successful_smp) == len(smp2lib):
            result['batchSchedule']['note'] = '分析完成'
        elif len(failed_smp):
            result['batchSchedule']['note'] = '分析异常'
            result['label'] = 'Error'
        log_text = f"{str(len(successful_smp))} 个样本分析完成\n" \
                   f"{str(len(failed_smp))} 个样本分析异常\n"
        result['batchSchedule']['finished'] = len(successful_smp)
        result['batchSchedule']['finished_samples'] = list(successful_smp)
        result['batchSchedule']['failed'] = len(failed_smp)
        result['batchSchedule']['failed_samples'] = list(failed_smp)
        result['batchSchedule']['processing'] = len(smp2lib) - len(successful_smp) - len(failed_smp)
        result['batchSchedule']['processing_samples'] = [s for s in smp2lib.keys() if s not in (successful_smp | failed_smp)]
        result['batchSchedule']['log'] = log_text
    else:
        # sample stat
        for smp in result['sampleSchedule'].keys():
            is_complete = True
            for ana in result['sampleSchedule'][smp].keys():
                if ana != '07.Stat' and result['sampleSchedule'][smp][ana]['note'] != '分析完成':
                    is_complete = False
            if is_complete:
                result['sampleSchedule'][smp]['07.Stat']['note'] = '分析中'
        # batch stat
        if len(failed_smp):
            result['batchSchedule']['note'] = '分析异常'
            result['batchSchedule']['failed'] = len(failed_smp)
            result['batchSchedule']['failed_samples'] = list(failed_smp)
            result['batchSchedule']['processing'] = len(smp2lib) - len(failed_smp)
            result['batchSchedule']['processing_samples'] = [s for s in smp2lib.keys() if s not in failed_smp]
        else:
            result['batchSchedule']['note'] = '分析中'
            result['batchSchedule']['processing'] = len(smp2lib)
            result['batchSchedule']['processing_samples'] = list(smp2lib.keys())

    return result
    # with open('haha.txt', 'w') as wfp:
    #     wfp.write('====================== 达安mNGS分析流程xxx批次分析日志 ======================\n\n')
    #     wfp.write(f'====================== 日志生成时间 {list(data.keys())[0]} ======================\n')
    #     wfp.write('\n')
    #     wfp.write('====================== 批次分析进度 ======================\n')
    #     wfp.write('批次样本总数\t分析完成样本数\t正在分析样本数\t分析异常样本数\n')
    #     wfp.write('\t'.join([str(result['batchSchedule']['total']),
    #                          str(result['batchSchedule']['finished']),
    #                          str(result['batchSchedule']['processing']),
    #                          str(result['batchSchedule']['failed'])]) + '\n')
    #     wfp.write('\n')
    #     wfp.write('====================== 样本分析进度 ======================\n')
    #     wfp.write('样本编号\t00.Merge\t01.Fastp\t02.RemoveHost\t03.Taxonomy-Classify\t04.Taxonomy-Coverage\t'
    #               '05.ResistantGene\t06.VF\t07.Stat\n')
    #     for smp in smp2lib.keys():
    #         note_lst = [smp]
    #         for ana in result['sampleSchedule'][smp].keys():
    #             note_lst.append(result['sampleSchedule'][smp][ana]['note'])
    #         wfp.write('\t'.join(note_lst) + '\n')
    #     wfp.write('\n')
    #     wfp.write('====================== 详细日志 ======================\n')
    #     for smp in smp2lib.keys():
    #         wfp.write(smp + ':\n')
    #         for ana in result['sampleSchedule'][smp].keys():
    #             wfp.write(f'  {ana}: \n')
    #             wfp.write(f"    note: {result['sampleSchedule'][smp][ana]['note']}\n")
    #             wfp.write(f"    label: {result['sampleSchedule'][smp][ana]['label']}\n")
    #             wfp.write(f"    log: \n{result['sampleSchedule'][smp][ana]['log']}")


@click.command(no_args_is_help=True)
@click.option('-s', '--samples', help='Sample list. <str>')
@click.option('-w', '--workdir', help='Where the PIP out dir is. <str>')
def main(samples, workdir):
    """ 检查各个步骤的分析状态，反馈分析进度; 每分钟更新一次日志直到by2.done生成 """
    while True:
        done_file = Path(workdir) / 'by2.done'
        if done_file.exists():
            break
        else:
            time_tag = time.strftime('%Y-%m-%d %H:%M:%S')
            res = workflow(samples, workdir, False)
            data = {time_tag: res}
            with open(Path(workdir) / 'shell/pipeline.log', 'a') as wfp, \
                    open(Path(workdir) / 'shell/pipeline.log.json', 'a') as wfp2:
                json.dump(data, wfp, ensure_ascii=False)
                wfp.write('\n')
                json.dump(data, wfp2, ensure_ascii=False, indent=4)
                wfp2.write('\n')
            time.sleep(60)

    time_tag = time.strftime('%Y-%m-%d %H:%M:%S')
    res = workflow(samples, workdir, True)
    data = {time_tag: res}
    with open(Path(workdir) / 'shell/pipeline.log', 'a') as wfp, \
            open(Path(workdir) / 'shell/pipeline.log.json', 'a') as wfp2:
        json.dump(data, wfp, ensure_ascii=False)
        wfp.write('\n')
        json.dump(data, wfp2, ensure_ascii=False, indent=4)
        wfp2.write('\n')


if __name__ == '__main__':
    main()