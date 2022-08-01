raw_data=/biodata/databases/autocnv_database/raw_data
version_data=/biodata/databases/autocnv_database/raw_data/2021-10-28
result_data=/biodata/databases/autocnv_database/data
mkdir -p $raw_data $version_data $result_data

#download raw data
wget -P $result_data https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
wget -P $result_data https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
wget -P $result_data https://ftp.ncbi.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz
wget -P $result_data https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget -P $result_data https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv
wget -P $result_data https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv
wget -P $result_data https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget -P $result_data https://decipher.sanger.ac.uk/files/downloads/HI_Predictions_Version3.bed.gz
wget -P $result_data https://datasetgnomad.blob.core.windows.net/dataset/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
wget -P $result_data https://datasetgnomad.blob.core.windows.net/dataset/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.bed.gz
wget -P $result_data ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/csv/genefamily_db_tables/family.csv

wget -P $raw_data https://www.omim.org/static/omim/data/mim2gene.txt
echo "gene_symbol" >$raw_data/gene-list-key-lte3.xlsx
grep gene $raw_data/mim2gene.txt |awk '{print $4}'|sort|uniq|sed '1d' >>$raw_data/gene-list-key-lte3.xlsx

wget -P $raw_data http://dgv.tcag.ca/dgv/docs/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.gff3

tabix -f $version_data/clinvar.vcf.gz

python liumm_create_database.py

#sort exon.sorted.bed
head -1 $result_data/exon.sorted.bed >temp.bed
bedtools sort -i $result_data/exon.sorted.bed >>temp.bed
mv temp.bed $result_data/exon.sorted.bed

head -1 manul_cnv_syndrome.txt >$result_data/cnv-syndrome-del.bed
cat manul_cnv_syndrome.txt |awk -F "\t" '$7~/Deletion/'|sort -k 1,1 -k 2,2n >>$result_data/cnv-syndrome-del.bed
head -1 manul_cnv_syndrome.txt >$result_data/cnv-syndrome-dup.bed
cat manul_cnv_syndrome.txt |awk -F "\t" '$7~/Duplication/'|sort -k 1,1 -k 2,2n >>$result_data/cnv-syndrome-dup.bed


#prepare bed files
ls $result_data/*bed $result_data/func-region.sorted|while read L
do
	sh bgzip_tabix.sh bed $L
done

#prepare vcf files
sed -i 's/AF_.\+=nan;//' $result_data/clinvar-pathogenic.sorted.vcf
ls $result_data/*vcf|while read L
do
	sh bgzip_tabix.sh vcf $L
done
