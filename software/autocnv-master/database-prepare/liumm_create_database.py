from multiprocessing import Pool
import pandas as pd
import numpy as np
import allel
import pysam
import json
import os
from datetime import datetime
pd.set_option('display.max_columns', None)
from urllib import request
from tqdm import tqdm
from autocnv import settings
from gtfparse import read_gtf

work_dir = settings.BASE_DIR
raw_data_dir = os.path.join(work_dir, 'raw_data')
result_data_dir = os.path.join(work_dir, 'data')
version_date = '2021-10-28'
version_date_dir = os.path.join(raw_data_dir, version_date)

# ori file
exon_ori_file = os.path.join(version_date_dir, 'hg19.refGene.gtf.gz')
cyto_band_ori_file = os.path.join(version_date_dir, 'cytoBand.txt.gz')
gene_info_ori_file = os.path.join(version_date_dir, 'Homo_sapiens.gene_info.gz')
ref_gene_ori_file = os.path.join(version_date_dir, 'refGene.txt.gz')
clingen_gene_ori_file = os.path.join(version_date_dir, 'ClinGen_gene_curation_list_GRCh37.tsv')
clingen_region_ori_file = os.path.join(version_date_dir, 'ClinGen_region_curation_list_GRCh37.tsv')
#clingen_gene_ori_file = os.path.join(version_date_dir, 'clingen_gene_hg19.tsv')
#clingen_region_ori_file = os.path.join(version_date_dir, 'clingen_region_hg19.tsv')
clinvar_ori_vcf_file = os.path.join(version_date_dir, 'clinvar.vcf.gz')
hi_pred_ori_file = os.path.join(version_date_dir, 'HI_Predictions_Version3.bed.gz')
gnomad_lof_ori_file = os.path.join(version_date_dir, 'gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
gnomad_control_ori_file = os.path.join(version_date_dir, 'gnomad_v2.1_sv.controls_only.sites.bed.gz')
# prep an external omim gene list
omim_gene_list_file = os.path.join(raw_data_dir, 'gene-list-key-lte3.xlsx')
dgv_ori_file = os.path.join(raw_data_dir, 'DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3')
# result file

gene_file = os.path.join(result_data_dir, 'gene.sorted.bed')
omim_gene_file = os.path.join(result_data_dir, 'omim-gene.sorted.bed')
clinvar_file = os.path.join(result_data_dir, 'clinvar-pathogenic.sorted.vcf')
decipher_gene_file = os.path.join(result_data_dir, 'decipher-gene.sorted.bed')
dgv_gain_file = os.path.join(result_data_dir, 'dgv-gain.sorted.bed')
dgv_loss_file = os.path.join(result_data_dir, 'dgv-loss.sorted.bed')
func_region_file = os.path.join(result_data_dir, 'func-region.sorted')
gnomad_del_file = os.path.join(result_data_dir, 'gnomad-del.sorted.bed')
gnomad_dup_file = os.path.join(result_data_dir, 'gnomad-dup.sorted.bed')
hi_cds_file = os.path.join(result_data_dir, 'hi-cds.sorted.bed')
hi_exon_file = os.path.join(result_data_dir, 'hi-exon.sorted.bed')
hi_gene_file = os.path.join(result_data_dir, 'hi-gene.sorted.bed')
hi_region_file = os.path.join(result_data_dir, 'hi-region.sorted.bed')
ts_gene_file = os.path.join(result_data_dir, 'ts-gene.sorted.bed')
ts_region_file = os.path.join(result_data_dir, 'ts-region.sorted.bed')
uhi_gene_file = os.path.join(result_data_dir, 'uhi-gene.sorted.bed')
uhi_region_file = os.path.join(result_data_dir, 'uhi-region.sorted.bed')
uts_gene_file = os.path.join(result_data_dir, 'uts-gene.sorted.bed')
uts_region_file = os.path.join(result_data_dir, 'uts-region.sorted.bed')
hgnc_gene_fam_file = os.path.join(version_date_dir, 'family.csv')


#1. create cyto-band.bed.gz
cyto_band_df = pd.read_csv(
    cyto_band_ori_file, sep='\t',
    names=['#chrom', 'start', 'end', 'name', 'gieStain']
)
cyto_band_df.to_csv(settings.CYTO_BAND_FILE.strip('.gz'), sep='\t', index=False)


#2. create gene.sorted.bed.gz
cols = [
    'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
    'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'ExonFrames']
refgene = pd.read_csv(ref_gene_ori_file, sep='\t', names=cols)
refgene = refgene[~refgene['chrom'].str.match(r'.*fix$')]
refgene['length'] = refgene['cdsEnd'] - refgene['cdsStart']
refgene = refgene.sort_values('length', ascending=False)
refgene = refgene.drop_duplicates('name2', keep='first')

gene_info = pd.read_csv(gene_info_ori_file, sep='\t')
refgene_info = refgene.merge(gene_info, left_on='name2', right_on='Symbol')
gene_fam = pd.read_csv(hgnc_gene_fam_file, sep=',')
gene_fam = gene_fam[gene_fam['typical_gene'].notna()].rename(columns={'id': 'fam'})
gene_fam_grp = gene_fam.groupby('typical_gene').agg({'fam': set})
gene = refgene_info.loc[
    refgene_info['type_of_gene'] == 'protein-coding',
    ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
].sort_values(['chrom', 'txStart', 'txEnd'])
gene.rename(columns={
    'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(gene_file, index=False, sep='\t')



#3. create omim-gene.sorted.bed.gz
raw_omim_df = pd.read_csv(omim_gene_list_file)
omim_gene = set(raw_omim_df['gene_symbol'].dropna())
omim_gene = refgene_info.loc[
    refgene_info['name2'].isin(omim_gene),
    ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
].sort_values(['chrom', 'txStart', 'txEnd'])
omim_gene.rename(columns={
    'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(omim_gene_file, index=False, sep='\t')


#4. 
curation_gene = pd.read_csv(clingen_gene_ori_file, sep='\t', dtype=str, skiprows=5)
hi_genes = set(
    curation_gene.loc[
        curation_gene['Haploinsufficiency Score'] == '3', '#Gene Symbol'
    ]
)
hi_gene = refgene_info.loc[
    refgene_info['name2'].isin(hi_genes),
    ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
].sort_values(['chrom', 'txStart', 'txEnd'])
hi_gene.rename(columns={
    'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(hi_gene_file, index=False, sep='\t')


#5.
hi_cds = refgene_info.loc[
    (refgene_info['name2'].isin(hi_genes)) & (refgene_info['length'] != 0),
    ['chrom', 'cdsStart', 'cdsEnd', 'GeneID', 'name2', 'name']
].sort_values(['chrom', 'cdsStart', 'cdsEnd'])
hi_cds.rename(columns={
    'chrom': '#chrom', 'cdsStart': 'start', 'cdsEnd': 'end',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(hi_cds_file, index=False, sep='\t')

#6.
hi_exons = refgene_info.loc[
    refgene_info['name2'].isin(hi_genes), ['chrom', 'exonStarts', 'exonEnds', 'GeneID', 'name2', 'name', 'strand']
].copy()
hi_exons['exonStarts'] = hi_exons['exonStarts'].str.replace(r',$', '')
hi_exons['exonEnds'] = hi_exons['exonEnds'].str.replace(r',$', '')
start = hi_exons['exonStarts'].str.split(',').apply(pd.Series).stack().reset_index()
start = start.rename(columns={'level_0': 'row', 0: 'start'})[['row', 'start']]
start['start'] = start['start'].astype(int)
end = hi_exons['exonEnds'].str.split(',').apply(pd.Series).stack().reset_index()
end = end.rename(columns={0: 'end'})['end'].astype(int)
position = start.join(end)
exon = position.merge(
    hi_exons[['chrom', 'GeneID', 'name2', 'name', 'strand']], how='left', left_on='row', right_index=True
)
exon = exon.sort_values(['chrom', 'start', 'end'])
exon['+'] = exon.groupby(['name2', 'name'])['start'].rank('first', ascending=True).astype(int)
exon['-'] = exon.groupby(['name2', 'name'])['start'].rank('first', ascending=False).astype(int)
exon['exon'] = pd.concat([exon.loc[exon['strand'] == '+', '+'], exon.loc[exon['strand'] == '-', '-']])
exon['last_exon'] = exon.groupby(['name2', 'name'])['exon'].transform('max') == exon['exon']
exon = exon[
    ['chrom', 'start', 'end', 'GeneID', 'name2', 'name', 'exon', 'last_exon']
].sort_values(['chrom', 'start', 'end'])
exon.rename(columns={
    'chrom': '#chrom',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(hi_exon_file, index=False, sep='\t')


#7.create clinvar-pathogenic.sorted.vcf
last_exon = exon[exon['last_exon']]
last_exon_region = last_exon['chrom'] + ':' + last_exon['start'].astype(str) + '-' + last_exon['end'].astype(str)
last_exon_region = last_exon_region.str.replace('chr', '')
need_fields = [
    'variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT',
    'variants/AF_ESP', 'variants/AF_EXAC', 'variants/AF_TGP', 'variants/CLNSIG'
]

with open(clinvar_file, 'w') as f:
    headers = allel.read_vcf_headers(clinvar_ori_vcf_file)
    f.write(''.join(headers.headers))

    def fetch_variants(region):
        fields, samples, headers, it = allel.iter_vcf_chunks(
            clinvar_ori_vcf_file, fields=need_fields, alt_number=1, region=region
        )
        for variants, *_ in it:
            esp_filter = np.isnan(variants['variants/AF_ESP'])
            esp_filter[~esp_filter] |= variants['variants/AF_ESP'][~esp_filter] < 0.01

            exac_filter = np.isnan(variants['variants/AF_EXAC'])
            exac_filter[~exac_filter] |= variants['variants/AF_EXAC'][~exac_filter] < 0.01

            tgp_filter = np.isnan(variants['variants/AF_TGP'])
            tgp_filter[~tgp_filter] |= variants['variants/AF_TGP'][~tgp_filter] < 0.01

            pathogenic_filter = np.isin(
                variants['variants/CLNSIG'], ['Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic']
            )

            af_filter = esp_filter & exac_filter & tgp_filter & pathogenic_filter

            filtered_variants = {k: v[af_filter] for k, v in variants.items()}

            filtered_variants['variants/CHROM'] = 'chr' + filtered_variants['variants/CHROM']

            return allel.normalize_callset(filtered_variants)

    with Pool(processes=32) as pool:
        variants = pool.map(fetch_variants, last_exon_region)

    for names, callset in tqdm(filter(lambda x: x is not None, variants)):
        allel.write_vcf_data(f, names, callset, None, {field: np.nan for field in need_fields})


#8.
uhi_genes = set(
    curation_gene.loc[curation_gene['Haploinsufficiency Score'] == '40', '#Gene Symbol']
)
uhi_gene = refgene_info.loc[
    refgene_info['name2'].isin(uhi_genes),
    ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
].sort_values(['chrom', 'txStart', 'txEnd'])
uhi_gene.rename(columns={
    'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(uhi_gene_file, index=False, sep='\t')


#9.
region = pd.read_csv(clingen_region_ori_file, sep='\t', skiprows=5, dtype=str)
position = region['Genomic Location'].str.extract(r'(?P<chrom>chr\w+)\s*:\s*(?P<start>\d+)\s*-\s*(?P<end>\d+)')
position['start'] = position['start'].astype(int)
position['end'] = position['end'].astype(int)
region_pos = region.merge(position, how='left', left_index=True, right_index=True)
func_region = region_pos.loc[
    (region_pos['Haploinsufficiency Score'].isin(['1', '2','3']))
     | (region_pos['Triplosensitivity Score'].isin(['1', '2','3'])),
     ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
].sort_values(['chrom', 'start', 'end'])
func_region.rename(columns={
    'chrom': '#chrom',
    '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
}).to_csv(func_region_file, sep='\t', index=False)

#10.
hi_region = region_pos.loc[
    region_pos['Haploinsufficiency Score'] == '3',
    ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
].sort_values(['chrom', 'start', 'end'])
hi_region['omim_genes'] = hi_region.apply(
    lambda row: ','.join(omim_gene.loc[
        (omim_gene['chrom'] == row['chrom'])
        & (omim_gene['txEnd'] >= row['start'])
        & (omim_gene['txStart'] <= row['end']),
        'name2'
    ]),
    axis=1
)
hi_region.rename(columns={
    'chrom': '#chrom',
    '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
}).to_csv(hi_region_file, sep='\t', index=False)

#11.
uhi_region = region_pos.loc[
    region_pos['Haploinsufficiency Score'] == '40',
    ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
].sort_values(['chrom', 'start', 'end'])
uhi_region['genes'] = uhi_region.apply(
    lambda row: ','.join(gene.loc[
        (gene['chrom'] == row['chrom'])
        & (gene['txEnd'] >= row['start'])
        & (gene['txStart'] <= row['end']),
        'name2'
    ]),
    axis=1
)
uhi_region.rename(columns={
    'chrom': '#chrom',
    '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
}).to_csv(uhi_region_file, sep='\t', index=False)

#12
decipher = pd.read_csv(hi_pred_ori_file, sep='\t',skiprows=1, header=None, usecols=[3,])
decipher = decipher[3].str.split('|', expand=True).rename(columns={0: 'symbol', 1: 'hi_score', 2: 'hi_index'})
decipher['hi_index'] = decipher['hi_index'].str.replace('%', '').astype(float)
decipher = decipher.merge(gene, left_on='symbol', right_on='name2')
gnomad = pd.read_csv(gnomad_lof_ori_file, sep='\t', index_col=0, compression='gzip')
decipher = decipher.join(gnomad['pLI'], on='name2')
decipher = decipher.join(gnomad['oe_lof_upper'], on='name2')
decipher = decipher.loc[
    (decipher['pLI'] >= 0.9) & (decipher['hi_index'] < 10) & (decipher['oe_lof_upper'] < 0.35),
    ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'pLI', 'hi_score']
].sort_values(['chrom', 'txStart', 'txEnd'])
decipher.rename(columns={
    'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(decipher_gene_file, sep='\t', index=False)


#13.
ts_genes = set(
    curation_gene.loc[
        curation_gene['Triplosensitivity Score'] == '3', '#Gene Symbol'
    ]
)
ts_gene = refgene_info.loc[
    refgene_info['name2'].isin(ts_genes),
    ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
].sort_values(['chrom', 'txStart', 'txEnd'])
ts_gene.rename(columns={
    'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(ts_gene_file, index=False, sep='\t')

#14.
uts_genes = set(
    curation_gene.loc[curation_gene['Triplosensitivity Score'] == '40', '#Gene Symbol']
)
uts_gene = refgene_info.loc[
    refgene_info['name2'].isin(uts_genes),
    ['chrom', 'txStart', 'txEnd', 'GeneID', 'name2', 'name', 'strand']
].sort_values(['chrom', 'txStart', 'txEnd'])
uts_gene.rename(columns={
    'chrom': '#chrom', 'txStart': 'start', 'txEnd': 'end',
    'GeneID': 'gene_id', 'name2': 'symbol', 'name': 'transcript'
}).to_csv(uts_gene_file, index=False, sep='\t')

#15.
ts_region = region_pos.loc[
    region_pos['Triplosensitivity Score'] == '3',
    ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
].sort_values(['chrom', 'start', 'end'])
ts_region['omim_genes'] = ts_region.apply(
    lambda row: ','.join(omim_gene.loc[
        (omim_gene['chrom'] == row['chrom'])
        & (omim_gene['txEnd'] >= row['start'])
        & (omim_gene['txStart'] <= row['end']),
        'name2'
    ]),
    axis=1
)
ts_region.rename(columns={
    'chrom': '#chrom',
    '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
}).to_csv(ts_region_file, sep='\t', index=False)

#16.
uts_region = region_pos.loc[
    region_pos['Triplosensitivity Score'] == '40',
    ['chrom', 'start', 'end', '#ISCA ID', 'ISCA Region Name']
].sort_values(['chrom', 'start', 'end'])
uts_region['genes'] = uts_region.apply(
    lambda row: ','.join(gene.loc[
        (gene['chrom'] == row['chrom'])
        & (gene['txEnd'] >= row['start'])
        & (gene['txStart'] <= row['end']),
        'name2'
    ]),
    axis=1
)
uts_region.rename(columns={
    'chrom': '#chrom',
    '#ISCA ID': 'isca_id', 'ISCA Region Name': 'name'
}).to_csv(uts_region_file, sep='\t', index=False)

#17.
dgv = pd.read_csv(
    dgv_ori_file,
    sep='\t', names=['chrom', 'info'], usecols=[0, 8]
).drop_duplicates('info')
info = dgv['info'].str.extract(
    r'ID=(?P<id>[^;]+).*variant_sub_type=(?P<type>[^;]+).*inner_start=(?P<start>[^;]+).*inner_end=(?P<end>[^;]+).*Frequency=(?P<freq>\S+?)%;.*num_unique_samples_tested=(?P<sample>[^;]+)'
).astype({'start': int, 'end': int, 'sample': int, 'freq': float})
dgv = dgv.merge(info, left_index=True, right_index=True)
dgv['af'] = dgv['freq'] / 100
dgv = dgv[dgv['sample'] >= 1000].sort_values(['chrom', 'start', 'end'])
def fetch_gene(gene, chrom, start, end):
    return ','.join(gene.loc[
        (gene['chrom'] == chrom) & (gene['txEnd'] >= start) & (gene['txStart'] <= end), 'name2'
    ])
with Pool(processes=7) as pool:
    dgv['genes'] = pool.starmap(fetch_gene, (
        (gene, row['chrom'], row['start'], row['end']) for _, row in dgv.iterrows()
    ), chunksize=70)
dgv.loc[
    dgv['type'] == 'Gain', ['chrom', 'start', 'end', 'id', 'genes', 'af']
].rename(columns={'chrom': '#chrom'}).to_csv(dgv_gain_file, sep='\t', index=False)


#18.
dgv.loc[
    dgv['type'] == 'Loss', ['chrom', 'start', 'end', 'id', 'genes', 'af']
].rename(columns={'chrom': '#chrom'}).to_csv(dgv_loss_file, sep='\t', index=False)


#19.
gnomad = pd.read_csv(
    gnomad_control_ori_file, sep='\t',  dtype=str,
    usecols=[0, 1, 2, 3, 4, 37, 38, 73, 74, 107, 108, 141, 142, 175, 176, 241]
)
gnomad = gnomad[
    (gnomad['FILTER'] == 'PASS') & gnomad['svtype'].isin(['DEL', 'DUP'])
]
gnomad = gnomad[
    (gnomad['N_BI_GENOS'].astype(int) >= 1000) |
    (gnomad['AFR_N_BI_GENOS'].astype(int) >= 1000) |
    (gnomad['AMR_N_BI_GENOS'].astype(int) >= 1000) |
    (gnomad['EAS_N_BI_GENOS'].astype(int) >= 1000) |
    (gnomad['EUR_N_BI_GENOS'].astype(int) >= 1000)
]
gnomad['#chrom'] = 'chr' + gnomad['#chrom']
gnomad['start'] = gnomad['start'].astype(int)
gnomad['end'] = gnomad['end'].astype(int)
with Pool(processes=30) as pool:
    gnomad['genes'] = pool.starmap(fetch_gene, (
        (gene, row['#chrom'], row['start'], row['end']) for _, row in gnomad.iterrows()
    ), chunksize=70)
gnomad['gene'] = gnomad.apply(lambda row: fetch_gene(gene, row['#chrom'], row['start'], row['end']), axis=1)
gnomad[
    gnomad['svtype'] == 'DEL'
][['#chrom', 'start', 'end', 'name', 'genes', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF']].rename(columns={
    'AF': 'af', 'AFR_AF': 'af_afr', 'AMR_AF': 'af_amr', 'EAS_AF': 'af_eas', 'EUR_AF': 'af_eur'
}).to_csv(gnomad_del_file, sep='\t', index=False)

#20.
gnomad.loc[
    gnomad['svtype'] == 'DUP',
    ['#chrom', 'start', 'end', 'name', 'genes', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF']
].rename(columns={
    'AF': 'af', 'AFR_AF': 'af_afr', 'AMR_AF': 'af_amr', 'EAS_AF': 'af_eas', 'EUR_AF': 'af_eur'
}).to_csv(gnomad_dup_file, sep='\t', index=False)

#21.
refGene_df = read_gtf(exon_ori_file)
refGene_df = refGene_df[refGene_df['feature'] == 'exon']
gene_db = pd.read_csv(settings.GENE_DATABASE, sep='\t')
merge_df = gene_db.merge(refGene_df, left_on='transcript', right_on='transcript_id')[[
    '#chrom', 'start_y', 'end_y', 'gene_id_x', 'symbol', 'exon_number', 'transcript'
]]
merge_df.columns = [x.strip('_x').strip('_y') for x in merge_df.columns]

merge_df.to_csv(settings.GENE_EXON_DATABASE.strip('.gz'), sep='\t', index=False)
