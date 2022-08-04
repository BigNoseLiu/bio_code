#!/usr/bin/env python
import gzip
import os
import sys
import click
import warnings
import pandas as pd
from pathlib import Path
Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg

warnings.filterwarnings('ignore')

############################ developer info ##############################
__AUTHOR__ = "Sujiawei"
__EMAIL__ = "su_jiawei163@163.com"
__VERSION__ = "V1.0"
__DATE__ = "2022-06-01"
##########################################################################

KINGDOM2TAG = {
    'Bacteria': 'Bacteria',
    'Eukaryota:Fungi': 'Eukaryota_Fungi',
    'Eukaryota:Parasite': 'Eukaryota_Parasite',
    'Eukaryota:Protozoa': 'Eukaryota_Protozoa',
    'SpecialPathogens': 'SpecialPathogens',
    'Viruses': 'Viruses'
}


def check_file(file):
    if not Path(file).exists():
        raise FileExistsError(f'{str(file)} not exists!')
    return


############################### Annotating by DB ###############################
def get_tax_info(db):
    """ get seq_taxid annotation from PIDB """
    infile1 = Path(db) / 'all.Taxonomy.txt'
    infile2 = Path(db) / 'GroupList.info'
    for file in infile1, infile2:
        check_file(file)

    tid2spe_tid = {}
    tid2kingdom = {}
    tid2name = {}
    gns_tid_set = set()
    grp_tid_set = set()
    spe_tid2grp_tid = {}

    df = pd.read_csv(infile1, sep='\t', dtype=str)
    for row in df.itertuples():
        tid2kingdom[row[3]] = row[5]
        tid2kingdom[row[4]] = row[5]
        tid2kingdom[row[17]] = row[5]
        tid2name[row[17]] = row[10]
        gns_tid_set.add(row[14])
        if row[5] == 'Viruses':
            tid2name[row[1]] = row[13]
            tid2name[row[3]] = row[13]
            tid2name[row[4]] = row[13]
            tid2spe_tid[row[3]] = row[1]
            tid2spe_tid[row[4]] = row[1]
        else:
            tid2name[row[1]] = row[11]
            tid2name[row[3]] = row[11]
            tid2name[row[4]] = row[11]
            tid2spe_tid[row[3]] = row[4]
            tid2spe_tid[row[4]] = row[4]

    df_grp = pd.read_csv(infile2, sep='\t', dtype=str)
    df_grp = df_grp[df_grp['Focus'] == '*']
    for row in df_grp.itertuples():
        grp_tid_set.add(row[1])
        spe_tid2grp_tid[row[3]] = row[1]
        tid2kingdom[row[1]] = row[9]
        tid2name[row[1]] = row[2]

    return {
        'tid2spe_tid': tid2spe_tid,
        'tid2kingdom': tid2kingdom,
        'tid2name': tid2name,
        'gns_tid_set': gns_tid_set,
        'grp_tid_set': grp_tid_set,
        'spe_tid2grp_tid': spe_tid2grp_tid
    }


@click.command(no_args_is_help=True)
@click.option('-m', '--mark', help='Enter a list of sequences id to validate. <str>')
@click.option('-k', '--kout', help='Kraken2 output file. <str>')
@click.option('--db', help='Path of PIDB <str>')
@click.option('--nt', help='Version of Nt Database <str> [20210103]', default="20210103")
@click.option('--fq', help='Clean data\'s fastq file, separated with comma, example: --fq fq1,fq2 <str>')
@click.option('-o', '--outdir', help='Results output dir <str>')
def main(mark, kout, db, nt, fq, outdir):
    """ Split sequence by tax """
    # check file
    check_file(mark)
    for file in fq.split(','):
        check_file(file)

    # get tax info from PIDB
    tax_info = get_tax_info(db)

    # get seqid to taxid
    seqid2taxid = {}
    with open(kout, 'r') as rfp:
        for _, line in enumerate(rfp):
            if not line.startswith('C'):
                continue
            ll = line.split('\t')
            if ll[2] not in ['0', '1', '2', '9606']:
                seqid2taxid[ll[1]] = ll[2]

    # deal with seq file
    data_seq = []
    for fq_file in fq.split(','):
        if fq_file.endswith('_1.fq.gz'):
            reads_tag = '1'
        elif fq_file.endswith('_2.fq.gz'):
            reads_tag = '2'
        else:
            reads_tag = '1'
        with gzip.open(fq_file, 'rb') as rfp:
            n = 0
            for line in rfp:
                n += 1
                if n % 4 == 1:
                    ll = str(line, 'utf-8').strip().split(' ')
                    seqid = ll[0].replace('@', '', 1)
                    if seqid in tax_info['grp_tid_set']:
                        taxid = seqid
                    else:
                        taxid = seqid2taxid.get(seqid)
                    item = [seqid, taxid]
                elif n % 4 == 2:
                    item.append(str(line, 'utf-8').strip())
                    item.append(reads_tag)
                    data_seq.append(item)
                    del item

    df = pd.DataFrame(data_seq, columns=['SeqID', 'Taxid', 'Seq', 'ReadsTag'])
    del data_seq

    seqid2taxid = {}
    with open(mark, 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.strip().split('\t')
            if ll:
                seqid2taxid[ll[1]] = ll[2]

    vld_fna = Path(outdir) / 'validated.fna'
    with open(vld_fna, 'w') as wfp:
        df_temp = df[df['SeqID'].isin(seqid2taxid)]
        # output for SE
        if len(fq.split(',')) == 1:
            for row in df_temp.itertuples():
                wfp.write(f'>{row[1]}|kraken:taxid|{row[2]}\n{row[3]}\n')
        # output for PE
        else:
            for row in df_temp.itertuples():
                wfp.write(f'>{row[1]}_{row[4]}|kraken:taxid|{row[2]}\n{row[3]}\n')

    # blast validate
    vld_m8 = Path(outdir) / 'validated.m8'
    cmd = f"{cfg.blastn} -query {str(vld_fna)} -db {cfg.BLASTDB}/{nt}/nt -num_threads 20 -out {str(vld_m8)} -outfmt " \
          f"\"6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp\" " \
          f"-evalue 0.05 -word_size 28 -reward 1 -penalty -2 -perc_identity 95 -qcov_hsp_perc 90"
    status = os.system(cmd)

    # deal with blast result
    seqid2score = {}
    with open(vld_m8, 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.strip().split('\t')
            temp_score = float(ll[2]) + float(ll[14])
            if seqid2score.get(ll[0], 0) < temp_score:
                seqid2score[ll[0]] = temp_score

    mark_seqid = set()
    with open(vld_m8, 'r') as rfp:
        for _, line in enumerate(rfp):
            ll = line.strip().split('\t')
            query_tid = ll[0].split('|')[-1]
            if ':' in ll[1]:
                subject_tid = ll[1].split(':')[0]
            elif '|' in ll[1]:
                subject_tid = ll[1].split('|')[0]
            else:
                continue
            if float(ll[2]) + float(ll[14]) == seqid2score.get(ll[0], '-'):
                if tax_info['tid2spe_tid'].get(query_tid, 'q') == tax_info['tid2spe_tid'].get(subject_tid, 't'):
                    mark_seqid.add(ll[0])
                if query_tid in tax_info['grp_tid_set']:
                    sbj_grp_tid = tax_info['spe_tid2grp_tid'].get(tax_info['tid2spe_tid'].get(subject_tid))
                    if query_tid == sbj_grp_tid:
                        mark_seqid.add(ll[0])

    vld_out = Path(outdir) / 'validated.out'
    with open(vld_out, 'w') as wfp, open(vld_fna, 'r') as rfp:
        for seqid in seqid2score.keys():
            tid = seqid.split('|')[-1]
            if tid in tax_info['grp_tid_set']:
                spe_name = tax_info['tid2name'].get(tid)
            else:
                spe_tid = tax_info['tid2spe_tid'].get(tid)
                spe_name = tax_info['tid2name'].get(spe_tid)
            if seqid in mark_seqid:
                if spe_name:
                    wfp.write(f'{seqid}\t{spe_name}\ttrue\n')
            else:
                if spe_name:
                    wfp.write(f'{seqid}\t{spe_name}\tfalse\n')
        # Fix some reads that missing verification results
        for _, line in enumerate(rfp):
            if line.startswith('>'):
                seqid = line.strip().replace('>', '', 1)
                tid = line.strip().split('|')[-1]
                if tid in tax_info['grp_tid_set']:
                    spe_name = tax_info['tid2name'].get(tid)
                else:
                    spe_tid = tax_info['tid2spe_tid'].get(tid)
                    spe_name = tax_info['tid2name'].get(spe_tid)
                if seqid not in seqid2score and spe_name:
                    wfp.write(f'{seqid}\t{spe_name}\tfalse\n')


if __name__ == '__main__':
    main()