#!/usr/bin/env python
import os
import re
import sys
import json
import click
import warnings
import pandas as pd
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

############################ developer info ##############################
__AUTHOR__  = "Sujiawei"
__EMAIL__   = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__    = "2022-05-06"
##########################################################################


warnings.filterwarnings('ignore')

KINGDOM2TAG = {
    'Bacteria': 'Bacteria',
    'Eukaryota:Fungi': 'Eukaryota_Fungi',
    'Eukaryota:Parasite': 'Eukaryota_Parasite',
    'Eukaryota:Protozoa': 'Eukaryota_Protozoa',
    'SpecialPathogens': 'SpecialPathogens',
    'Viruses': 'Viruses'
}


def check_file(file):
    if not file.exists():
        raise FileExistsError(f'{str(file)} not exists!')
    return


def get_tax_info(db, sample_type):
    """ 获取数据库中物种的注释信息 """
    tax_file = Path(db) / 'all.Taxonomy.txt'
    check_file(tax_file)

    tax_infos = defaultdict(defaultdict)
    type2index = {"RT": 19, "Plasma": 20, "CSF": 21, "Others": 22}

    df = pd.read_csv(tax_file, sep='\t', dtype=str)
    for row in df.itertuples():
        if row[5] == 'Viruses':
            tax_infos['name2Kingdom'][row[13]] = row[5]
            tax_infos['name2Genus'][row[13]] = row[10]
            tax_infos['name2Type'][row[13]] = row[2]
            tax_infos['name2ChineseName'][row[13]] = row[18]
            tax_infos['name2SampleSite'][row[13]] = row[type2index.get(sample_type)]
            tax_infos['name2Info'][row[13]] = row[24]
            tax_infos['name2Medicine'][row[13]] = row[25]
            tax_infos['name2Reference'][row[13]] = row[26]
            tax_infos['name2Type_Judge'][row[13]] = row[27]
            tax_infos['taxid2name'][row[3]] = row[13]
            tax_infos['name2len'][row[13]] = row[14]
        else:
            tax_infos['name2Kingdom'][row[11]] = row[5]
            tax_infos['name2Genus'][row[11]] = row[10]
            tax_infos['name2Type'][row[11]] = row[2]
            tax_infos['name2ChineseName'][row[11]] = row[18]
            tax_infos['name2SampleSite'][row[11]] = row[type2index.get(sample_type)]
            tax_infos['name2Info'][row[11]] = row[24]
            tax_infos['name2Medicine'][row[11]] = row[25]
            tax_infos['name2Reference'][row[11]] = row[26]
            tax_infos['name2Type_Judge'][row[11]] = row[27]
            tax_infos['taxid2name'][row[3]] = row[11]
            tax_infos['taxid2name'][row[4]] = row[11]
            tax_infos['name2len'][row[11]] = row[15]
        tax_infos['taxid2kingdom'][row[3]] = row[5]
        tax_infos['taxid2kingdom'][row[4]] = row[5]
        tax_infos['genus2GenusCN'][row[10]] = row[23]
        if row[16] == '*':
            tax_infos['patho'][row[3]] = 1
            tax_infos['patho'][row[4]] = 1

    # Group infos for target groups
    grp_file = Path(db) / 'GroupList.info'
    check_file(grp_file)
    df_grp = pd.read_csv(grp_file, sep='\t', dtype=str)
    df_grp = df_grp[df_grp['Focus'] == '*']
    for row in df_grp.itertuples():
        tax_infos['name2Kingdom'][row[2]] = row[9]
        tax_infos['name2Genus'][row[2]] = row[7]
        tax_infos['name2Info'][row[2]] = row[11]
        tax_infos['name2Medicine'][row[2]] = row[12]
        tax_infos['name2Reference'][row[2]] = row[13]
        tax_infos['name2Type'][row[2]] = row[14]
        tax_infos['name2Type_Judge'][row[2]] = row[15]
        tax_infos['name2ChineseName'][row[2]] = row[5]
        tax_infos['taxid2name'][row[1]] = row[2]
        tax_infos['spe2group'][row[4]] = row[2]
        tax_infos['patho'][row[1]] = 1

    # group taxid set
    grp_tid_set = set(df_grp['GroupTaxid'])

    # group to species set
    grp2spe_set = defaultdict(set)
    for row in df_grp[['GroupName', 'SpeciesName']].itertuples():
        grp2spe_set[row[1]].add(row[2])

    # all tax
    all_tid_set = set(df['Taxid']) | set(df['SpeciesTaxid']) | grp_tid_set

    # pathogen's names set
    df_temp = df[df['Focus'] == '*']
    patho_names = set(df_temp['Species']) | set(df_temp['OrganismName']) | \
                  set(df_temp['SubtypeName']) | set(df_grp['GroupName'])

    tax_infos['grp2spe_set'] = grp2spe_set
    tax_infos['grp_tid_set'] = grp_tid_set
    tax_infos['all_tid_set'] = all_tid_set
    tax_infos['patho_names'] = patho_names

    return tax_infos


def get_validate_result(infile):
    """ 获取blast验证结果 """
    nodelete_species = set()
    nodelete_groups = set()
    spe2check = {}
    with open(infile, 'r') as rfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                if ll[2] == 'true':
                    nodelete_species.add(ll[1])
                    spe2check[ll[1]] = 'True'
                    if 'group' in ll[1] or 'complex' in ll[1]:
                        nodelete_groups.add(ll[1])

    delete_species = set()
    deleted_groups = set()
    with open(infile, 'r') as rfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                if ll[1] not in nodelete_species:
                    delete_species.add(ll[1])
                    if 'group' in ll[1] or 'complex' in ll[1]:
                        deleted_groups.add(ll[1])
    validated_res = {
        'nodelete_species': nodelete_species,
        'nodelete_groups': nodelete_groups,
        'spe2check': spe2check,
        'delete_species': delete_species,
        'deleted_groups': deleted_groups
    }
    return validated_res


def get_score(highconf):
    """ 获取taxid对应的reads最高kmer得分 """
    species2score = {}
    with open(highconf, 'r') as rfp:
        for line in rfp.readlines():
            line = line.strip()
            if line:
                ll = line.split('\t')
                species2score[ll[0]] = ll[1]
    return species2score


def stat_summary(name, mpa, kreport, outdir):
    """ Summary """
    out_tax = Path(outdir) / 'Taxonomy_Summary.txt'

    stat = {"Viruses": 0, "Bacteria": 0, "Archaea": 0, "Fungi": 0,
            "Metazoa": 0, "Eukaryota": 0, "Human": 0, "Unclassified": 0,
            "Classified": 0, "Total_Reads": 0, "Unclassified_Rate": 0,
            "Classified_Rate": 0, "Metazoa_Parasite": 0, "Protozoa": 0}

    with open(mpa, 'r') as rfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                if ll[0] == 'd__Viruses':
                    stat['Viruses'] = int(ll[1])
                if ll[0] == 'd__Bacteria':
                    stat['Bacteria'] = int(ll[1])
                if ll[0] == 'd__Archaea':
                    stat['Archaea'] = int(ll[1])
                if ll[0] == 'd__Eukaryota|k__Fungi':
                    stat['Fungi'] = int(ll[1])
                if ll[0] == 'd__Eukaryota|k__Metazoa':
                    stat['Metazoa'] = int(ll[1])
                if ll[0] == 'd__Eukaryota':
                    stat['Eukaryota'] = int(ll[1])
                if re.search('\|g__Homo$', ll[0]):
                    stat['Human'] = int(ll[1])

    spe2reads_cnt = {}

    with open(kreport, 'r') as rfp, open(out_tax, 'w') as wfp:
        wfp.write("Sample\tTotal_Reads\tUnclassified_Reads_Number\tUnclassified_Rate\tClassified_Reads_Number\t"
                  "Classified_Rate\tViruses\tBacteria\tArchaea\tFungi\tProtozoa\tMetazoa_Parasite\tHuman\n")
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                species = ll[5].strip()
                spe2reads_cnt[species] = ll[1]
                if species == 'unclassified':
                    stat['Unclassified'] = int(ll[1])
                if species == 'root':
                    stat['Classified'] = int(ll[1])
        stat['Total_Reads'] = stat['Unclassified'] + stat['Classified']
        if stat['Total_Reads'] == 0:             # err
            stat['Unclassified_Rate'] = 0
            stat['Classified_Rate'] = 0
        else:
            stat['Unclassified_Rate'] = format(stat['Unclassified'] / stat['Total_Reads'] * 100, ".2f")
            stat['Classified_Rate'] = format(stat['Classified'] / stat['Total_Reads'] * 100, ".2f")
        stat['Metazoa_Parasite'] = stat['Metazoa'] - stat['Human']
        stat['Protozoa'] = stat['Eukaryota'] - stat['Fungi'] - stat['Metazoa']
        wfp.write(f"{name}\t{stat['Total_Reads']}\t{stat['Unclassified']}\t{stat['Unclassified_Rate']}\t"
                  f"{stat['Classified']}\t{stat['Classified_Rate']}")
        wfp.write(f"\t{stat['Viruses']}\t{stat['Bacteria']}\t{stat['Archaea']}\t{stat['Fungi']}\t"
                  f"{stat['Protozoa']}\t{stat['Metazoa_Parasite']}\t{stat['Human']}\n")

    return spe2reads_cnt


def custom_dict():
    return {
        'is_subtype': False,
        'preSpecies': '',
        'all_mapped': 0,   # unique and multiple mapped reads count
        'uni_mapped': 0,   # unique mapped reads count
    }


def zero():
    return 0


def stat_tax(kreport, tax_infos, validate_res, spe2score, outdir):
    """ 统计 kreport """
    out_taxid = Path(outdir) / 'taxid.list'

    sum_dict = {'sum_reads': 0, 'uniform': 0}
    spe_level2rc = defaultdict()
    grp_level2rc = defaultdict()
    virus2stat = defaultdict(custom_dict)

    preClass = ''
    preSpecies = ''
    with open(kreport, 'r') as rfp, open(out_taxid, 'w') as wfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if not line:
                continue
            ll = line.split('\t')
            kingdom = tax_infos['taxid2kingdom'].get(ll[4])
            species = tax_infos['taxid2name'].get(ll[4])

            ll[1] = int(ll[1])
            ll[2] = int(ll[2])

            if ll[3] == 'S':
                preClass = kingdom
                preSpecies = species

            # stat for target groups
            if ll[4] in tax_infos['grp_tid_set']:
                if spe2score.get(species) or species not in validate_res['deleted_groups']:
                    sum_dict['sum_reads'] += ll[1]
                    grp_level2rc[species] = ll[2]

            # stat for species level
            if ll[3] == 'S' and ll[4] != '9606':
                if spe2score.get(species) or species not in validate_res['delete_species']:
                    sum_dict['sum_reads'] += ll[1]
                    if ll[4] in tax_infos['all_tid_set']:
                        wfp.write(f'\n{preClass}\t{preSpecies}\t{ll[4]}')
                        if kingdom == 'Viruses':
                            virus2stat[species]['all_mapped'] += ll[1]
                            virus2stat[species]['uni_mapped'] += ll[2]
                        else:
                            spe_level2rc[species] = ll[1]

            # stat for subtype level of viruses
            if re.match('^S\d', ll[3]) and (ll[4] in tax_infos['all_tid_set']):
                if spe2score.get(species) or species not in validate_res['delete_species']:
                    wfp.write(f',{ll[4]}')
                if preClass == 'Viruses':
                    virus2stat[species]['all_mapped'] += ll[1]
                    virus2stat[species]['uni_mapped'] += ll[2]
                    if species != preSpecies:
                        virus2stat[species]['is_subtype'] = True
                        virus2stat[species]['preSpecies'] = preSpecies

    for virus, stat in virus2stat.items():
        if stat['is_subtype']:
            if spe2score.get(virus) or virus not in validate_res['delete_species']:
                if virus not in spe_level2rc:
                    spe_level2rc[virus] = stat['all_mapped']
                else:
                    spe_level2rc[virus] += stat['all_mapped']
        else:
            if stat['uni_mapped'] != 0:
                if virus not in spe_level2rc:
                    spe_level2rc[virus] = stat['uni_mapped']
                else:
                    spe_level2rc[virus] += stat['all_mapped']

    # 统计检出物种中，属(Genus)最大reads数
    gns2spe_rc_lst = defaultdict(list)
    for spe, spe_rc in spe_level2rc.items():
        genus = tax_infos['name2Genus'].get(spe)
        if genus and genus != '-':
            gns2spe_rc_lst[genus].append(spe_rc)
    gns2gns_max = {}
    for gns, lst in gns2spe_rc_lst.items():
        gns2gns_max[gns] = max(lst)

    # species RPM
    spe2rpm = defaultdict()
    for spe in spe_level2rc.keys():
        if sum_dict['sum_reads'] == 0:         # err
            spe2rpm[spe] = 0
        else:
            spe2rpm[spe] = round(spe_level2rc[spe] / sum_dict['sum_reads'] * 10 ** 6, 2)

    # stat for group reads count
    group2rc = defaultdict(zero)
    for grp, spe_set in tax_infos['grp2spe_set'].items():
        group2rc[grp] += grp_level2rc.get(grp, 0)
        for spe in spe_set:
            if spe2score.get(spe) or validate_res['spe2check'].get(spe) == 'True':
                group2rc[grp] += spe_level2rc.get(spe, 0)

    # group RPM
    grp2rpm = {}
    for grp, grp_rc in group2rc.items():
        if sum_dict['sum_reads'] == 0:
            grp2rpm[grp] = 0
        else:
            grp2rpm[grp] = round(grp_rc / sum_dict['sum_reads'] * 10 ** 6, 2)

    # stat for group genomes length
    group2len = {}
    for grp, spe_set in tax_infos['grp2spe_set'].items():
        if grp not in group2len:
            group2len[grp] = 0
        for spe in spe_set:
            spe_len = tax_infos['name2len'].get(spe)
            if spe_len and re.match('^[1-9]', spe_len):
                group2len[grp] += int(spe_len)

    # prepare for RH-RPM, species
    spe2uniform = {}
    for spe in spe2rpm.keys():
        spe_len = tax_infos['name2len'].get(spe)
        if spe_len and re.match('^[1-9]', spe_len):
            value = spe_level2rc[spe] / int(spe_len)
        else:
            print(f'Warning: genome sequence length of {spe} not exists.')
            value = spe_level2rc[spe] / 1
        spe2uniform[spe] = value
        sum_dict['uniform'] += value
    # prepare for RH-RPM, group
    grp2uniform = {}
    for grp in grp2rpm.keys():
        grp_rc = grp_level2rc.get(grp, 0)
        grp_len = group2len.get(grp, 0)
        if grp_len:
            value = grp_rc / grp_len
        else:
            print(f'Warning: genome sequence length of {grp} not exists.')
            value = grp_rc / 1
        grp2uniform[grp] = value
        sum_dict['uniform'] += value

    # stat for RH-RPM
    # species
    spe2tpm = {}
    for spe, uniform in spe2uniform.items():
        if sum_dict['uniform'] == 0:
            spe2tpm[spe] = 0
        else:
            spe2tpm[spe] = round(uniform / sum_dict['uniform'] * 10 ** 6, 2)
    # group
    grp2tpm = {}
    for grp, uniform in grp2uniform.items():
        if sum_dict['uniform'] == 0:
            grp2tpm[grp] = 0
        else:
            grp2tpm[grp] = round(uniform / sum_dict['uniform'] * 10 ** 6, 2)

    return {
        'spe2rpm': spe2rpm,
        'spe2tpm': spe2tpm,
        'grp2rpm': grp2rpm,
        'grp2tpm': grp2tpm,
        'genus2genus_max': gns2gns_max,
        'spe_level2rc': spe_level2rc,
        'group2rc': group2rc,
    }


def modify_kingdom(df_data):
    """ 云康定制化需求, 修改特定物种的Kingdom信息便于LIMS展示 """
    # 展示属中所有的物种
    target_gns_set = {'Brucella', 'Legionella', 'Mycobacterium', 'Nocardia', }
    # 当检出指定物种时, 展示指定物种对应属的所有物种
    target_spe2gns = {
        'Yersinia pestis': 'Yersinia',
        'Vibrio cholerae': 'Vibrio',
        'Vibrio vulnificus': 'Vibrio',
        'Bacillus anthracis': 'Bacillus',
        'Listeria monocytogenes': 'Listeria',
        'Cryptococcus neoformans': 'Cryptococcus',
        'Pneumocystis jirovecii': 'Pneumocystis',
        'Talaromyces marneffei': 'Talaromyces'
    }

    df_sub = df_data[df_data['SpeciesSN'].isin(target_spe2gns)]
    if not df_sub.empty:
        for gns in set(df_sub['GenusSN']):
            target_gns_set.add(gns)

    new_kd_list = []
    old_kd_list = []
    for row in df_data[['Kingdom', 'GenusSN', 'SpeciesSN']].itertuples():
        old_kd_list.append(row[1])
        if row[2] in target_gns_set:
            new_kd_list.append('SpecialPathogens')
        else:
            new_kd_list.append(row[1])
    df_data['Kingdom'] = new_kd_list
    df_data['rawKingdom'] = old_kd_list

    return df_data


def sort_filter(df_data, target_pathogen, top=10):
    """ 按界进行归类，再按属reads、种reads分别由高到低进行排序; 根据致病菌分级列表进行筛选, 每个属最多保留前top个物种结果 """
    """ 统计属排位 """
    df_tmp = pd.DataFrame(columns=df_data.columns)
    for kd in set(df_data['Kingdom']):
        df_sub = df_data[df_data['Kingdom'] == kd]
        gns_rc2rank = {}
        for i, v in enumerate(sorted(set(df_sub['GenusReads']), reverse=True), 1):
            gns_rc2rank[v] = i
        df_sub['GenusRank'] = [gns_rc2rank.get(rc) for rc in df_sub['GenusReads']]
        df_tmp = pd.concat([df_tmp, df_sub])
    df_data = df_tmp

    """ 物种排序及过滤 """
    df_ret = pd.DataFrame(columns=df_data.columns)
    for kd in set(df_data['Kingdom']):
        df_kd = df_data[df_data['Kingdom'] == kd]
        df_kd = df_kd.astype({'GenusReads': int, 'SpeciesReads': int})
        df_sort = df_kd.sort_values(by=['GenusReads', 'SpeciesReads'], ascending=False)
        df_filter = pd.DataFrame(columns=df_data.columns)
        for gns in set(df_sort['GenusSN']):
            df_gns = df_sort[df_sort['GenusSN'] == gns]
            # is pathogens
            df_is_ptg = df_gns[df_gns['SpeciesSN'].isin(target_pathogen)]
            df_grp = df_gns[df_gns['GroupSN'] == df_gns['SpeciesSN']]
            df_filter = pd.concat([df_filter, df_is_ptg, df_grp])
            df_filter = df_filter.drop_duplicates(subset=['SpeciesSN'])
            if len(df_is_ptg) >= 10:
                continue
            # is not pathogens
            df_not_ptg = df_gns[~df_gns['SpeciesSN'].isin(target_pathogen)]
            if df_not_ptg.empty:
                continue
            df_not_ptg = df_not_ptg.sort_values(by=['SpeciesReads'], ascending=False)
            df_filter = pd.concat([df_filter, df_not_ptg.iloc[:top-len(df_is_ptg)]])
        df_filter = df_filter.sort_values(by=['GenusReads', 'SpeciesReads'], ascending=False)
        df_ret = pd.concat([df_ret, df_filter])
    return df_ret


def sort(df_data):
    """ 按界进行归类，再按属reads、种reads分别由高到低进行排序 """
    df_ret = pd.DataFrame(columns=df_data.columns)
    for kd in set(df_data['Kingdom']):
        df_kd = df_data[df_data['Kingdom'] == kd]
        df_kd = df_kd.astype({'GenusReads': int, 'SpeciesReads': int})
        df_sort = df_kd.sort_values(by=['GenusReads', 'SpeciesReads'], ascending=False)
        for gns in set(df_sort['GenusSN']):
            df_gns = df_sort[df_sort['GenusSN'] == gns]
            df_ret = pd.concat([df_ret, df_gns])
    return df_ret


TXT_HEADER = "Kingdom\tGenusSN\tGenusCN\tGenusReads\tGenusReads(Normalization)\tGenusRelativeAbundance(%)\t" \
             "GenusRank\tGroupSN\tGroupCN\tGroupReads\tGroupReads(Normalization)\tSpeciesSN\tSpeciesCN\t" \
             "SpeciesReads\tSpeciesReads(Normalization)\tSpeciesRelativeAbundance(%)\tSpeciesRank\t" \
             "GMRN\tTPM\tRH-RPM\tVerification\tCoverage\tBinCoverage\tFocus\tLabel\tType\tDatabaseCode\t" \
             "Freq/[Min-25%-Median-75%-Max](Historically)\tFoldChange(Historically)\tNCCount\tFoldChange(NC)\t" \
             "Background_level\n"


def write_data(wfp, spe_tid, spe, spe_lv_rc, anno_info, anno_mark, check_tag, score):
    wfp.write(f"{spe_tid}\t{spe}\t{spe_lv_rc}\t-\t{anno_info}")
    if check_tag:
        wfp.write(f"\t{check_tag}\t{anno_mark}\n")
    else:
        if score:
            wfp.write(f"\t{score}\t{anno_mark}\n")
        else:
            wfp.write(f"\t-\t{anno_mark}\n")
    return wfp


@click.command(no_args_is_help=True)
@click.option('--report', help="kraken's report file")
@click.option('--mpa', help="kraken's mpa-report file")
@click.option('--db', help="path of the database to annotate")
@click.option('--name', help="sample name")
@click.option('--validate', help="Mark2Fa.pl's output")
@click.option('--highconf', help="High confidence kmer score")
@click.option('--type', help="RT|Plasma|CSF|[Others]")
@click.option('--outdir', help="output path. [./]")
@click.option('--version', help="print version information.", is_flag=True)
def main(report, mpa, db, name, validate, highconf, type, outdir, version):
    """ Stat Taxonomy results """
    if version:
        print(f"\tScript Version:   {__VERSION__}")
        print(f"\tLast Modified At: {__DATE__}")
        return

    if not Path(outdir).exists():
        os.system(f'mkdir -p {outdir}')

    for file in report, mpa, validate, highconf:
        # 检查输入文件是否存在
        check_file(Path(file))

    tax_infos = get_tax_info(db, type)

    validate_res = get_validate_result(validate)

    spe2reads_cnt = stat_summary(name, mpa, report, outdir)

    spe2score = get_score(highconf)

    tax_stat_res = stat_tax(report, tax_infos, validate_res, spe2score, outdir)

    # output tax stat results #
    file_total = Path(outdir).joinpath('Total_Detail.txt')

    with open(file_total, 'w') as out_total:
        out_total.write(TXT_HEADER)
        for spe, rpm in tax_stat_res['spe2rpm'].items():
            tpm = str(tax_stat_res['spe2tpm'].get(spe))
            rpm = str(tax_stat_res['spe2rpm'].get(spe))
            superkingdom = tax_infos['name2Kingdom'].get(spe)
            genus = tax_infos['name2Genus'].get(spe)
            genus_cn = tax_infos['genus2GenusCN'].get(genus, '-')
            spe_cname = tax_infos['name2ChineseName'].get(spe, '-')
            label = tax_infos['name2Type_Judge'].get(spe, '-')
            spe_type = tax_infos['name2Type'].get(spe)

            genus_rc = str(spe2reads_cnt.get(genus, 0))
            genus_max_rc = str(tax_stat_res['genus2genus_max'].get(genus, 0))
            group = tax_infos['spe2group'].get(spe, '-')
            group_cn = tax_infos['name2ChineseName'].get(group, '-')
            group_rc = str(tax_stat_res['group2rc'].get(group, 0))
            spe_lv_rc = str(tax_stat_res['spe_level2rc'].get(spe))

            target_tag = tax_infos['name2SampleSite'].get(spe)
            score = spe2score.get(spe)
            check_tag = validate_res['spe2check'].get(spe)
            if not check_tag and not score:
                continue

            if superkingdom:
                out_total.write(superkingdom)
                if genus == 'NA' or not genus:
                    out_total.write('\t-\t-\t0\t0\t0\t0')
                else:
                    out_total.write(f"\t{genus}\t{genus_cn}\t{genus_rc}\t0\t0\t0")
                if group:
                    out_total.write(f"\t{group}\t{group_cn}\t{group_rc}\t0")
                else:
                    out_total.write("\t-\t-\t0\t0")
                out_total.write(f"\t{spe}\t{spe_cname}\t{spe_lv_rc}\t0\t0\t0\t{genus_max_rc}\t{tpm}\t{rpm}")
                if check_tag:
                    out_total.write(f"\t{check_tag}\t0\t0\t{target_tag}\t{label}\t{spe_type}\t-\t0\t0\t0\t0\t-\n")
                else:
                    if score:
                        out_total.write(f"\t{score}\t0\t0\t{target_tag}\t{label}\t{spe_type}\t-\t0\t0\t0\t0\t-\n")
                    # else:
                    #     out_total.write(f"\tFalse\t0\t0\t{target_tag}\t{label}\t{spe_type}\t-\t0\t0\t0\t0\t-\n")

        # 输出复合群相关的信息
        for grp, grp_rc in tax_stat_res['group2rc'].items():
            grp_rc = tax_stat_res['group2rc'].get(grp, 0)
            if not grp_rc:
                continue
            else:
                grp_rc = str(grp_rc)
            skd = tax_infos['name2Kingdom'].get(grp, '-')
            gns = tax_infos['name2Genus'].get(grp, '-')
            gns_cn = tax_infos['genus2GenusCN'].get(gns, '-')
            gns_rc = str(spe2reads_cnt.get(gns, 0))
            gns_max_rc = str(tax_stat_res['genus2genus_max'].get(gns, 0))
            grp_cn = tax_infos['name2ChineseName'].get(grp, '-')
            grp_mark = tax_infos['name2Type'].get(grp, '-')
            grp_label = tax_infos['name2Type_Judge'].get(grp, '-')
            grp_tpm = str(tax_stat_res['grp2tpm'].get(grp, 0))
            grp_rpm = str(tax_stat_res['grp2rpm'].get(grp, 0))
            if spe2score.get(grp):
                verification = spe2score.get(grp)
            else:
                s_score_lst = []
                for s in tax_infos['grp2spe_set'].get(grp):
                    if spe2score.get(s):
                        s_score_lst.append(float(spe2score.get(s)))
                if s_score_lst:
                    verification = str(max(s_score_lst))
                else:
                    if validate_res['spe2check'].get(grp):
                        verification = validate_res['spe2check'].get(grp)
                    else:
                        verification = '-'
                        for s in tax_infos['grp2spe_set'].get(grp):
                            if validate_res['spe2check'].get(s):
                                verification = validate_res['spe2check'].get(s)
            # output for out_total, out_pathogen
            txt_row = [skd, gns, gns_cn, gns_rc, '0', '0', '0',
                       grp, grp_cn, grp_rc, '0',
                       grp, grp_cn, grp_rc, '0', '0', '0', gns_max_rc,
                       grp_tpm, grp_rpm, verification, '0\t0\t+', grp_label, grp_mark, '-\t0\t0\t0\t0\t-\n']

            out_total.write('\t'.join(txt_row))

    # stat reads count
    tmp_df = pd.read_csv(file_total, sep='\t', dtype=str)
    tmp_df = tmp_df.fillna('-')
    tmp_df = tmp_df.astype({'GenusReads': int, 'GroupReads': int, 'SpeciesReads': int, 'GMRN': int})
    # stat for GenusRelativeAbundance
    gns2rc = {row[1]: row[2] for row in tmp_df[['GenusSN', 'GenusReads']].itertuples()}
    sum_gns_rc = sum(gns2rc.values())
    gns_abundance = []
    for gns_rc in tmp_df['GenusReads']:
        if not sum_gns_rc:
            gns_abundance.append('0.000')
            continue
        g_abu = gns_rc / sum_gns_rc * 100
        if g_abu == 0:
            gns_abundance.append('0.000')
        elif g_abu < 0.001:
            gns_abundance.append('{:.3e}'.format(g_abu))
        else:
            gns_abundance.append('{:.3f}'.format(g_abu))
    tmp_df['GenusRelativeAbundance(%)'] = gns_abundance

    # stat for Genus Rank
    gns_rc2rank = {}
    for i, gns_rc in enumerate(sorted(set(gns2rc.values()), reverse=True), 1):
        gns_rc2rank[gns_rc] = i
    tmp_df['GenusRank'] = [gns_rc2rank.get(i) for i in tmp_df['GenusReads']]

    # stat for SpeciesRelativeAbundance
    sum_spe_rc = sum(tmp_df['SpeciesReads'])
    spe_abundance = []
    for spe_rc in tmp_df['SpeciesReads']:
        s_abu = spe_rc / sum_spe_rc * 100
        if s_abu == 0:
            spe_abundance.append('0.000')
        elif s_abu < 0.001:
            spe_abundance.append('{:.3e}'.format(s_abu))
        else:
            spe_abundance.append('{:.3f}'.format(s_abu))
    tmp_df['SpeciesRelativeAbundance(%)'] = spe_abundance

    # sort raw data
    tmp_df = sort(tmp_df)

    # output raw data
    tmp_df.to_csv(file_total.with_suffix('.raw.txt'), index=False, sep='\t')
    tmp_df.to_excel(file_total.with_suffix('.raw.xlsx'), index=False)

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

    # Modify species's kingdom
    tmp_df = modify_kingdom(tmp_df)

    # Sort and filter data
    tmp_df = sort_filter(tmp_df, target_pathogen)

    # output filtered data
    tmp_df.to_excel(file_total.with_suffix('.xlsx'), index=False)
    tmp_df.to_csv(file_total, index=False, sep='\t')

    file_pathogen = Path(outdir) / 'Pathogeny_Detail.txt'
    df_pathogen = tmp_df[tmp_df['SpeciesSN'].isin(tax_infos['patho_names'])]
    df_pathogen.to_csv(file_pathogen, sep='\t', index=False)
    df_pathogen.to_excel(file_pathogen.with_suffix('.xlsx'), index=False)


if __name__ == '__main__':
    main()