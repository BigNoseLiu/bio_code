#!/usr/bin/env python
import json
import os
import re
import sys
import click
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


def check_file(file):
    if not file.exists():
        raise FileExistsError(f'{str(file)} not exists!')
    return


# get annotation from db
def get_species_info(db):
    """ 获取物种注释信息 """
    infile = Path(db) / 'all.latin2chinese.INFO.txt'
    check_file(infile)
    df = pd.read_csv(infile, sep='\t', dtype=str)
    df2 = df.copy()
    df = df.drop_duplicates(subset=['SpeciesEN'], keep='last')
    df.index = df['SpeciesEN']
    df = df.fillna('None')
    spe2anno = df.to_dict('index')

    df2 = df2.drop_duplicates(subset=['Taxid'], keep='last')
    df2.index = df2['Taxid']
    df2 = df2.fillna('None')
    taxid2anno = df2.to_dict('index')
    return spe2anno, taxid2anno


def get_tax_info(db, sample_type):
    """ 获取数据库中物种分类的相关信息 """
    tax_infos = defaultdict(defaultdict)
    type2index = {"RT": 15, "Plasma": 16, "CSF": 17, "Others": 18}

    for file in ['all.Taxonomy.txt', 'all.Taxonomy.other.txt']:
        file_obj = Path(db) / file
        check_file(file_obj)
        fp = file_obj.open('r')
        next(fp)
        for _, line in enumerate(fp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                tax_infos['spe2class'][ll[8]] = ll[2]
                tax_infos['spe2class'][ll[9]] = ll[2]
                tax_infos['spe2spe_tid'][ll[8]] = ll[1]
                tax_infos['strain2seq_tid'][ll[9]] = ll[0]
                tax_infos['strain2len'][ll[9]] = ll[10]
                tax_infos['spe2len'][ll[8]] = ll[11]
                # tax_infos['spe2CN'][ll[8]] = ll[14]
                # tax_infos['spe2CN'][ll[9]] = ll[14]
                tax_infos['spe2genus'][ll[8]] = ll[7]
                tax_infos['spe2genus'][ll[9]] = ll[7]
                tax_infos['taxid2spe'][ll[1]] = ll[8]
                tax_infos['spe_tid2class'][ll[1]] = ll[2]
                if ll[2] == 'Viruses':
                    tax_infos['taxid2spe'][ll[0]] = ll[9]
                    tax_infos['spe2spe_tid'][ll[9]] = ll[1]
                    tax_infos['speciesLv'][ll[8]] = 1
                    tax_infos['seq_tid2spe_tid'][ll[0]] = ll[1]
                    tax_infos['spe_tid2speLv'][ll[1]] = ll[8]
                else:
                    tax_infos['taxid2spe'][ll[0]] = ll[8]
                tax_infos['spe2target'][ll[8]] = ll[type2index[sample_type]]
                tax_infos['spe2target'][ll[9]] = ll[type2index[sample_type]]
                if ll[12] == '*':
                    tax_infos['patho'][ll[8]] = 1
                    tax_infos['patho'][ll[9]] = 1
        fp.close()
    return tax_infos


def get_group_info(db):
    """ 获取复合群注释信息 """
    infile = Path(db) / 'GroupList.info'
    check_file(infile)
    df = pd.read_csv(infile, sep='\t', dtype=str)
    df = df.drop_duplicates(subset=['GroupName'])
    df = df.fillna('None')
    df.index = df['GroupName']
    group2infos = df.to_dict('index')

    infile2 = Path(db) / 'SpeciesGroup.list'
    check_file(infile2)
    spe2group = {}
    with open(infile2, 'r') as rfp:
        next(rfp)
        for _, line in enumerate(rfp):
            ll = line.strip().split('\t')
            if ll:
                spe2group[ll[3]] = ll[1]

    return group2infos, spe2group


def get_seqtid2spe2tid(db):
    infile = Path(db) / 'seqTaxid2speciesTaxid'
    check_file(infile)
    seqtid2spetid = defaultdict()
    with open(infile, 'r') as rfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                seqtid2spetid[ll[0]] = ll[1]
    return seqtid2spetid


def get_validate_result(infile, spe2group):
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
                    if ll[1] in spe2group:
                        nodelete_groups.add(spe2group[ll[1]])

    delete_species = set()
    deleted_groups = set()
    with open(infile, 'r') as rfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                if ll[1] not in nodelete_species:
                    delete_species.add(ll[1])
                    if ll[1] in spe2group and spe2group[ll[1]] not in nodelete_groups:
                        deleted_groups.add(spe2group[ll[1]])
    validated_res = {
        'nodelete_species': nodelete_species,
        'nodelete_groups': nodelete_groups,
        'spe2check': spe2check,
        'delete_species': delete_species,
        'deleted_groups': deleted_groups
    }
    return validated_res


def correct_group(kreport, spe2reads_cnt, delete_species, spe2group, group2infos):
    """ 根据validated.out的统计结果矫正复合群的检出reads"""
    # group2minus = {}
    group2reads_cnt = {}
    with open(kreport, 'r') as rfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if line:
                ll = line.split('\t')
                # species = ll[5].strip()
                sci_name = ll[5].strip()
                if sci_name in group2infos:
                    group2reads_cnt[sci_name] = int(ll[2])
                # if species in delete_species and species in spe2group:
                #     group = spe2group[species]
                #     if group not in group2minus:
                #         group2minus[group] = int(ll[1])
                #     else:
                #         group2minus[group] += int(ll[1])

    # for species in delete_species:
    #     if species in spe2group:
    #         # 去除blast验证为假的物种的键值
    #         del spe2group[species]

    # 从复合群reads数中减去验证为假的物种检出reads数
    # group2reads_cnt = {}
    # for group in group2infos.keys():
    #     if group in spe2reads_cnt:
    #         if group in group2minus:
    #             group2reads_cnt[group] = int(spe2reads_cnt[group]) - group2minus[group]
    #         else:
    #             group2reads_cnt[group] = spe2reads_cnt[group]
    # 过滤reads count = 0 的复合群
    # group4del = []
    # for group, rc in group2reads_cnt.items():
    #     if rc == 0:
    #         group4del.append(group)
    # for grp in group4del:
    #     del group2reads_cnt[grp]

    return group2reads_cnt, spe2group


def get_score(highconf, tax_infos):
    """ 获取taxid对应的reads最高kmer得分 """
    species2score = {}
    with open(highconf, 'r') as rfp:
        for line in rfp.readlines():
            line = line.strip()
            if line:
                ll = line.split('\t')
                if ll[0] in tax_infos['taxid2spe']:
                    species2score[tax_infos['taxid2spe'][ll[0]]] = ll[1]
                if ll[1] in tax_infos['taxid2spe']:
                    species2score[tax_infos['taxid2spe'][ll[0]]] = ll[1]
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
        'common': 0,
        'uniqTotal': 0,
        'types': defaultdict(),
        'S1_spe': set()
    }


def stat_tax(kreport, tax_infos, validate_res, spe2score, outdir):
    """ 统计 kreport """
    out_taxid = Path(outdir) / 'taxid.list'

    sum_dict ={'S': 0, 'uniform': 0}
    spe_level2rc = defaultdict()
    spe_level2stat = defaultdict(custom_dict)
    # 标记只输出S1水平的病毒
    virus2showS1 = defaultdict()
    genus2genus_max = defaultdict()
    spe1_2_spe = {}
    s1_lv = ''
    s1_to_sn = defaultdict(set)
    preClass = ''
    preSpecies = ''
    with open(kreport, 'r') as rfp, open(out_taxid, 'w') as wfp:
        for _, line in enumerate(rfp):
            line = line.strip()
            if not line:
                continue
            ll = line.split('\t')
            species = ll[5].strip()
            superkingdom = tax_infos['spe2class'].get(species)
            ll[1] = int(ll[1])
            ll[2] = int(ll[2])
            #
            if ll[3] == 'S':
                preClass = superkingdom
                preSpecies = species

            if ll[3] == 'S' and ll[4] != '9606':
                if spe2score.get(species) or species not in validate_res['delete_species']:
                    sum_dict['S'] += ll[1]
                    if superkingdom == 'Viruses':
                        virus2showS1[species] = defaultdict()
                        spe_level2stat[species]['common'] = ll[2]
                    else:
                        spe_level2rc[species] = ll[1]
                    if species in tax_infos['patho']:
                        wfp.write(f'\n{preClass}\t{preSpecies}\t{ll[4]}')

            if re.match('^S\d', ll[3]) and (preSpecies in tax_infos['patho']) \
                    and (preSpecies not in validate_res['delete_species']):
                wfp.write(f',{ll[4]}')

            if re.match('^S\d', ll[3]):
                if preClass == 'Viruses':
                    if spe2score.get(species) or preSpecies not in validate_res['delete_species']:
                        spe_level2stat[preSpecies]['uniqTotal'] += ll[2]
                        spe_level2stat[preSpecies]['types'][species] = ll[1]
                        if ll[3] == 'S1':
                            virus2showS1[preSpecies][species] = True
                            s1_lv = species
                            spe1_2_spe[species] = preSpecies
                            spe_level2stat[preSpecies]['S1_spe'].add(species)
                        else:
                            s1_to_sn[s1_lv].add(species)
            # 统计检出的非病毒物种中，属最大reads数
            if species in tax_infos['spe2genus'] and tax_infos['spe2class'][species] != 'Viruses':
                genus = tax_infos['spe2genus'][species]
                if genus not in genus2genus_max:
                    genus2genus_max[genus] = ll[1]
                elif genus2genus_max[genus] < ll[1]:
                    genus2genus_max[genus] = ll[1]

    # delete viruses that validated result is false
    spe_lv4del = set()
    for spe_lv in spe_level2stat.keys():
        if not spe2score.get(spe_lv) and spe_lv in validate_res['delete_species']:
            spe_lv4del.add(spe_lv)
        else:
            if len(spe_level2stat[spe_lv]['S1_spe']) > 0:
                # 过滤blast验证为假的亚种（S1水平）
                strain4del = set()
                for s1_spe in spe_level2stat[spe_lv]['S1_spe']:
                    if not spe2score.get(s1_spe) and s1_spe in validate_res['delete_species']:
                        strain4del.add(s1_spe)
                        del virus2showS1[spe_lv][s1_spe]
                for strain in strain4del:
                    if strain in spe_level2stat[spe_lv]['types']:
                        del spe_level2stat[spe_lv]['types'][strain]
                # 当所有S1水平物种验证为假且没有reads直接分类到S水平上，过滤该物种
                del_tag = True
                if spe_level2stat[spe_lv]['common'] > 0:
                    del_tag = False
                for strain in spe_level2stat[spe_lv]['S1_spe']:
                    if strain in spe_level2stat[spe_lv]['types']:
                        del_tag = False
                if del_tag:
                    spe_lv4del.add(spe_lv)
    for spe_lv in spe_lv4del:
        if spe_lv in spe_level2stat:
            del spe_level2stat[spe_lv]

    strain2species = defaultdict()
    for spe_lv in spe_level2stat.keys():
        if len(spe_level2stat[spe_lv]['types']) == 0:
            spe_level2rc[spe_lv] = spe_level2stat[spe_lv]['common']
            if spe_lv in tax_infos['spe2genus']:
                genus = tax_infos['spe2genus'][spe_lv]
                if (genus not in genus2genus_max) or (genus2genus_max[genus] < spe_level2rc[spe_lv]):
                    genus2genus_max[genus] = spe_level2rc[spe_lv]
        else:
            #############    病毒统计结果后续仍需优化    ##################
            for strain in spe_level2stat[spe_lv]['types']:
                if strain not in spe_level2stat[spe_lv]['S1_spe']:
                    continue
                strain_rc = spe_level2stat[spe_lv]['types'][strain]    #
                spe_uniqtotal = spe_level2stat[spe_lv]['uniqTotal']    #
                spe_common = spe_level2stat[spe_lv]['common']          #
                spe_level2rc[strain] = int(strain_rc / spe_uniqtotal * spe_common) + strain_rc
                strain2species[strain] = spe_lv
                if strain in tax_infos['spe2genus']:
                    genus = tax_infos['spe2genus'][strain]
                    if (genus not in genus2genus_max) or (genus2genus_max[genus] < spe_level2rc[strain]):
                          genus2genus_max[genus] = spe_level2rc[strain]

    spe2uniform = defaultdict()
    spe2RPM = defaultdict()
    for spe in spe_level2rc.keys():
        """ 以下代码块需要在更新all.Taxonomy.txt表中序列长度为空的项后再次进行修正 """
        ##########################################################################################
        value = 0
        if spe in strain2species:
            if spe in tax_infos['strain2len']:
                if re.match('^[1-9]', tax_infos['strain2len'][spe]):
                    value = spe_level2rc[spe] / float(tax_infos['strain2len'][spe])
                else:
                    value = spe_level2rc[spe] / 1
            else:
                if re.match('^[1-9]', tax_infos['spe2len'][strain2species[spe]]):
                    value = spe_level2rc[spe] / float(tax_infos['spe2len'][strain2species[spe]])
                else:
                    value = spe_level2rc[spe] / 1
        else:
            if spe in tax_infos['spe2len']:
                if re.match('^[1-9]', tax_infos['spe2len'][spe]):
                    value = spe_level2rc[spe] / float(tax_infos['spe2len'][spe])
                else:
                    value = spe_level2rc[spe] / 1
        ##########################################################################################
        if value and (spe2score.get(spe) or spe not in validate_res['delete_species']):
            spe2uniform[spe] = value
            spe2RPM[spe] = round(spe_level2rc[spe] / sum_dict['S'] * 1000000, 2)
            sum_dict['uniform'] += value

    return {
        'sum': sum_dict,
        'spe2uniform': spe2uniform,
        'spe2RPM': spe2RPM,
        'genus2genus_max': genus2genus_max,
        'spe_level2rc': spe_level2rc,
        'strain2species': strain2species,
        'S1_to_Sn': s1_to_sn,
        'virus2showS1': virus2showS1
    }


TXT_HEADER = "Kingdom\tGenusSN\tGenusCN\tGenusReads\tGenusReads(Normalization)\tGenusRelativeAbundance(%)\t" \
             "GenusRank\tGroupSN\tGroupCN\tGroupReads\tGroupReads(Normalization)\tSpeciesSN\tSpeciesCN\t" \
             "SpeciesReads\tSpeciesReads(Normalization)\tSpeciesRelativeAbundance(%)\tSpeciesRank\t" \
             "GMRN\tTPM\tRH-RPM\tVerification\tCoverage\tBinCoverage\tFocus\tLabel\tType\tDatabaseCode\t" \
             "Freq/[Min-25%-Median-75%-Max](Historically)\tFoldChange(Historically)\tNCCount\tFoldChange(NC)\t" \
             "Background_level\n"
# XLS_HEADER = "tax_id\tspecies\treads_count\tcoverage\tspecies_cn\tgenus_cn\tnoun\t" \
#              "medicine\tncbi\tblast_maxscore\tmark\ttype_judge\n"


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

    spe2anno, taxid2anno = get_species_info(db)

    # seqtid2spetid = get_seqtid2spe2tid(db)

    group2infos, spe2group = get_group_info(db)

    tax_infos = get_tax_info(db, type)

    validate_res = get_validate_result(validate, spe2group)
    #
    spe2reads_cnt = stat_summary(name, mpa, report, outdir)

    group2reads_cnt, spe2group = correct_group(report,
                                               spe2reads_cnt,
                                               validate_res['delete_species'],
                                               spe2group,
                                               group2infos)

    spe2score = get_score(highconf, tax_infos)

    tax_stat_res = stat_tax(report, tax_infos, validate_res, spe2score, outdir)

    # output tax stat results #
    file_total = Path(outdir).joinpath('Total_Detail.txt')
    file_pathogen = Path(outdir).joinpath('Pathogeny_Detail.txt')
    # file_bacteria = Path(outdir).joinpath('bacteria.xls')
    # file_fungi = Path(outdir).joinpath('fungi.xls')
    # file_ly = Path(outdir).joinpath('LY.xls')
    # file_protozoa = Path(outdir).joinpath('protozoa.xls')
    # file_viral = Path(outdir).joinpath('viral.xls')

    with open(file_total, 'w') as out_total, open(file_pathogen, 'w') as out_pathogen:
            # open(file_bacteria, 'w') as out_bacteria, open(file_fungi, 'w') as out_fungi, \
            # open(file_ly, 'w') as out_ly, open(file_protozoa, 'w') as out_protozoa, \
            # open(file_viral, 'w') as out_viral:
        out_total.write(TXT_HEADER)
        out_pathogen.write(TXT_HEADER)
        # out_bacteria.write(XLS_HEADER)
        # out_fungi.write(XLS_HEADER)
        # out_ly.write(XLS_HEADER)
        # out_protozoa.write(XLS_HEADER)
        # out_viral.write(XLS_HEADER)

        sort = defaultdict(dict)
        group_has_spe = set()
        for spe, rpm in tax_stat_res['spe2RPM'].items():
            tpm = format(tax_stat_res['spe2uniform'].get(spe) / tax_stat_res['sum']['uniform'] * 1000000, ".2f")
            rpm = str(tax_stat_res['spe2RPM'].get(spe))
            superkingdom = tax_infos['spe2class'].get(spe)
            genus = tax_infos['spe2genus'].get(spe)
            if spe in spe2anno:
                genus_cn = spe2anno[spe]['GenusCN']
                spe_cname = spe2anno[spe]["SpeciesCN"]
                label = spe2anno[spe]['Type_Judge']
                spe_type = spe2anno[spe]['[G+, G-]']
            else:
                genus_cn = '-'
                spe_cname = '-'
                label = '-'
                spe_type = '-'
            genus_rc = str(spe2reads_cnt.get(genus))
            genus_max_rc = str(tax_stat_res['genus2genus_max'].get(genus))
            group = spe2group.get(spe)
            if group:
                group_has_spe.add(group)
                group_cn = group2infos[group]['GroupChineseName']
            else:
                group_cn = '-'
            group_rc = str(group2reads_cnt.get(group))
            spe_lv_rc = str(tax_stat_res['spe_level2rc'].get(spe))

            spe_tid = tax_infos['spe2spe_tid'].get(spe)
            target_tag = tax_infos['spe2target'].get(spe)
            score = spe2score.get(spe)
            check_tag = validate_res['spe2check'].get(spe)
            if superkingdom == 'Viruses':
                if not score:
                    score = spe2score.get(tax_stat_res['strain2species'].get(spe))
                    if not score:
                        # 若S1水平病毒对应的Sn水平的病毒中有一个的blast验证为true,则check_tag为真,否则为假/None
                        if spe in tax_stat_res['S1_to_Sn']:
                            for sn in tax_stat_res['S1_to_Sn'][spe]:
                                check_tag = validate_res['spe2check'].get(sn)
                                if check_tag:
                                    break
            if not score and not check_tag:
                continue
            # if spe in spe2anno:
            #     anno_info = '\t'.join([spe2anno[spe]["SpeciesCN"], spe2anno[spe]["GenusCN"], spe2anno[spe]["Info"],
            #                            spe2anno[spe]["Medicine"], spe2anno[spe]["Reference"]])
            #     anno_mark = '\t'.join([spe2anno[spe]["[G+, G-]"], spe2anno[spe]["Type_Judge"]])
            # elif spe_tid in taxid2anno:
            #     anno_info = '\t'.join([taxid2anno[spe_tid]["SpeciesCN"], taxid2anno[spe_tid]["GenusCN"], taxid2anno[spe_tid]["Info"],
            #                            taxid2anno[spe_tid]["Medicine"], taxid2anno[spe_tid]["Reference"]])
            #     anno_mark = '\t'.join([taxid2anno[spe_tid]["[G+, G-]"], taxid2anno[spe_tid]["Type_Judge"]])
            # else:
            #     anno_info = '-\t-\t-\t-\t-'
            #     anno_mark = '-\t-'
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
                    else:
                        out_total.write(f"\tFalse\t0\t0\t{target_tag}\t{label}\t{spe_type}\t-\t0\t0\t0\t0\t-\n")

                if superkingdom != 'Viruses':
                    if not re.search(' sp\.| genomosp\.| genosp\.', spe):
                        out_pathogen.write(superkingdom)
                        if genus and genus != 'NA':
                            out_pathogen.write(f"\t{genus}\t{genus_cn}\t{genus_rc}\t0\t0\t0")
                        else:
                            out_pathogen.write('\t-\t-\t0\t0\t0\t0')
                        if group:
                            out_pathogen.write(f"\t{group}\t{group_cn}\t{group_rc}\t0")
                        else:
                            out_pathogen.write("\t-\t-\t0\t0")
                        out_pathogen.write(f"\t{spe}\t{spe_cname}\t{spe_lv_rc}\t0\t0\t0\t{genus_max_rc}\t{tpm}\t{rpm}")
                        if check_tag:
                            out_pathogen.write(f"\t{check_tag}\t0\t0\t{target_tag}\t{label}\t{spe_type}\t-\t0\t0\t0\t0\t-\n")
                        else:
                            if score:
                                out_pathogen.write(f"\t{score}\t0\t0\t{target_tag}\t{label}\t{spe_type}\t-\t0\t0\t0\t0\t-\n")
                            # else:
                            #     out_pathogen.write(f"\t-\t0\t-\t{target_tag}\t-\t-\t-\t0\t0\t0\t0\t-\n")
                    # if superkingdom == 'Bacteria':
                    #     out_bacteria = write_data(out_bacteria, spe_tid, spe, spe_lv_rc,
                    #                               anno_info, anno_mark, check_tag, score)
                    # elif superkingdom == 'Eukaryota:Fungi':
                    #     out_fungi = write_data(out_fungi, spe_tid, spe, spe_lv_rc,
                    #                            anno_info, anno_mark, check_tag, score)
                    # elif superkingdom == 'Eukaryota:Protozoa' or superkingdom == 'Eukaryota:Parasite':
                    #     out_protozoa = write_data(out_protozoa, spe_tid, spe, spe_lv_rc,
                    #                               anno_info, anno_mark, check_tag, score)
                    # else:
                    #     out_ly = write_data(out_ly, spe_tid, spe, spe_lv_rc,
                    #                         anno_info, anno_mark, check_tag, score)
                else:
                    if not re.search(' sp\.| genomosp\.| genosp\.', spe):
                        # 如果病毒有S1水平，则输出S1的结果；否则输出S水平的结果
                        strain2bool = tax_stat_res['virus2showS1'].get(spe)
                        if strain2bool and len(strain2bool) and spe in tax_infos['speciesLv']:
                            continue
                        out_pathogen.write(superkingdom)
                        if genus and genus != 'NA':
                            out_pathogen.write(f"\t{genus}\t{genus_cn}\t{genus_rc}\t0\t0\t0")
                        else:
                            out_pathogen.write('\t-\t-\t0\t0\t0\t0')
                        if group:
                            out_pathogen.write(f"\t{group}\t{genus_cn}\t{group_rc}\t0")
                        else:
                            out_pathogen.write("\t-\t-\t0\t0")
                        out_pathogen.write(f"\t{spe}\t{spe_cname}\t{spe_lv_rc}\t0\t0\t0\t{genus_max_rc}\t{tpm}\t{rpm}")
                        # strain_tid = str(tax_infos['strain2seq_tid'].get(spe))
                        if check_tag:
                            out_pathogen.write(f"\t{check_tag}\t0\t0\t{target_tag}\t{label}\t{spe_type}\t-\t0\t0\t0\t0\t-\n")
                            # out_viral.write(f"{strain_tid}\t{spe}\t{spe_lv_rc}\t-\t{anno_info}\t{check_tag}\t{anno_mark}\n")
                        else:
                            if score:
                                out_pathogen.write(f"\t{score}\t0\t0\t{target_tag}\t{label}\t{spe_type}\t-\t0\t0\t0\t0\t-\n")
                                # out_viral.write(f"{strain_tid}\t{spe}\t{spe_lv_rc}\t-\t{anno_info}\t{score}\t{anno_mark}\n")
                            # else:
                            #     out_pathogen.write(f"\tFalse\t{target_tag}\n")
                            #     out_viral.write(f"{strain_tid}\t{spe}\t{spe_lv_rc}\t-\t{anno_info}\tFalse\t{anno_mark}\n")
                species = tax_stat_res['strain2species'].get(spe)
                if species:
                    sort[superkingdom][species] = 1
                    sort[superkingdom][species] = len(sort[superkingdom].keys())
                else:
                    sort[superkingdom][spe] = 1
                    sort[superkingdom][spe] = len(sort[superkingdom].keys())

        # 输出复合群相关的信息
        for grp in group2reads_cnt.keys():
            if grp in group_has_spe:
                continue
            if group2infos[grp]['Focus'] == '*':
                gns = group2infos[grp]['Genus']
                gns_cn = group2infos[grp]['GenusCN']
                skd = group2infos[grp]['Superkingdom']
                grp_rc = str(group2reads_cnt.get(grp, 0))
                if grp_rc == '0':
                    continue
                gns_rc = str(spe2reads_cnt.get(gns, 0))
                gns_max_rc = str(tax_stat_res['genus2genus_max'].get(gns, 0))
                grp_cn = group2infos[grp]['GroupChineseName']
                # gns_cn = group2infos[grp]['GenusCN']
                # grp_mark = group2infos[grp]['Mark']
                # output for out_total, out_pathogen
                txt_row = [skd, gns, gns_cn, gns_rc, '0', '0', '0',
                           grp, grp_cn, grp_rc, '0',
                           grp, grp_cn, grp_rc, '0', '0', '0', gns_max_rc,
                           '0\t0\t-\t0\t0\t+', label, spe_type, '-\t0\t0\t0\t0\t-\n']
                out_total.write('\t'.join(txt_row))
                out_pathogen.write('\t'.join(txt_row))
                # output for out_bacteria, out_fungi, out_ly, out_protozoa, out_viral
                # grp_tid = int(group2infos[grp].get('GroupTaxid', 0))
                # xls_row = [str(grp_tid), grp, grp_rc, '-', grp_cn, gns_cn, group2infos[grp]['Info'], '-',
                #            group2infos[grp]['Reference'], '-', grp_mark, group2infos[grp].get('Type_Judge', '-')]
                # if skd == 'Viruses':
                #     out_viral.write('\t'.join(xls_row))
                # else:
                #     if skd == 'Bacteria':
                #         out_bacteria.write('\t'.join(xls_row) + '\n')
                #     elif skd == 'Eukaryota:Fungi':
                #         out_fungi.write('\t'.join(xls_row) + '\n')
                #     elif skd == 'Eukaryota:Protozoa' or skd == 'Eukaryota:Protozoa':
                #         out_protozoa.write('\t'.join(xls_row) + '\n')
                #     else:
                #         out_ly.write('\t'.join(xls_row) + '\n')

    # stat reads count
    for file in file_total, file_pathogen:
        tmp_df = pd.read_csv(file, sep='\t', dtype=str)
        tmp_df = tmp_df.fillna('-')
        tmp_df = tmp_df.astype({'GenusReads': int, 'GroupReads': int, 'SpeciesReads': int,
                                'GMRN': int, 'TPM': float, 'RH-RPM': float})
        # stat for GenusRelativeAbundance
        gns2rc = {row[1]: row[2] for row in tmp_df[['GenusSN', 'GenusReads']].itertuples()}
        sum_gns_rc = sum(gns2rc.values())
        gns_abundance = []
        for gns_rc in tmp_df['GenusReads']:
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
        tmp_df['GenusRank'] = tmp_df.apply(lambda x: gns_rc2rank.get(x['GenusReads']), axis=1)

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

        # output data
        tmp_df.to_excel(file.with_suffix('.xlsx'), index=False)
        tmp_df.to_csv(file.with_suffix('.tmp'), index=False, sep='\t')
        os.system(f"mv {str(file.with_suffix('.tmp'))} {str(file)}")

    taxid_r = Path(outdir).joinpath('taxid.list')
    taxid_w = Path(outdir).joinpath('taxid.list.tmp')
    with open(taxid_r, 'r') as rfp, open(taxid_w, 'w') as wfp:
        for i, line in enumerate(rfp, 1):
            line = line.strip()
            if line:
                ll = line.split('\t')
                if sort[ll[0]].get(ll[1]):
                    wfp.write(line + '\t' + str(sort[ll[0]].get(ll[1])) + '\n')
    os.system(f"mv {str(taxid_w)} {str(taxid_r)}")


if __name__ == '__main__':
    main()