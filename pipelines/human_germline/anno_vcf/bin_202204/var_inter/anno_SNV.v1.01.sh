#1. input files
bin_dir=$(dirname $0)
out_dir=$1

#2. system configs
docker_name=liumingming1988/biodocker
docker_cmd="docker run --rm -v /:/mnt -v $bin_dir/../../../../../:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u $ug_id $docker_name bash -c"

shift
annovar=/biodata/git_code/software/annovar/bin/table_annovar.pl
humandb=/mnt/data/databases/annovar_humandb/
factera=/biodata/git_code/software/factera/factera.pl
factera_exon=/biodata/git_code/software/factera/exons.no_chr.bed
ref=/biodata/databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta
ref_2bit=/biodata/databases/gatk_bundle/b37/human_g1k_v37_decoy.2bit








#&J	A202208020001_merge_all"
$docker_cmd "mkdir -p $out_dir"
$docker_cmd "perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/gatk4_germline_SNVIndel/merge_vcf.v1.01.pl  /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/germline_anno/cohort_test/cohort_test.final.leftalign.vcf /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/germline_anno/test/test.final.leftalign.vcf /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/germline_anno/test_chrm/test_chrm.final.leftalign.vcf >/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.vcf"
$docker_cmd "cat /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/exomiser/exome.yml|perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/exomiser/create_yml.v2.01.pl   liumm_ped_file_path:/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.fix.ped liumm_vcf_file_in:/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.vcf liumm_vcf_file_path:/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.exomiser_fix.vcf liumm_output_prefix:/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/exomiser.merge_all liumm_hpo_ids:HP:0000407,HP:0000365 >/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/exomiser.merge_all.yml"
$docker_cmd "java -Xms2g -Xmx4g -jar /biodata/git_code/software/exomiser/exomiser-cli-13.0.0/exomiser-cli-13.0.0.jar --analysis /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/exomiser.merge_all.yml --spring.config.location=/biodata/git_code/software/exomiser/exomiser-cli-13.0.0/application.properties"
$docker_cmd "java -jar /biodata/git_code/software/snpeff/snpEff5_0/snpEff.jar -v GRCh37.p13.RefSeq -s /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/snpeff.merge_all.stat.html /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.vcf >/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/snpeff.merge_all.vcf"
$docker_cmd "cat /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/snpeff.merge_all.vcf|perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/snpeff/fix_snpeff_Eff.v1.pl|perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/snpeff/select_standard_hgvs_from_snpeff_Eff.v2.01.pl|perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/gatk4_germline_SNVIndel/tag_nearbyOverlap_variant.v1.01.pl -ref /biodata/databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta >/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.vcf"
$docker_cmd "perl /biodata/git_code/software/annovar/bin/table_annovar.pl -buildver hg19 -protocol clinvar_20220416,clinvar_samepos_20220416,clinvar_flank5_20220416,HGMD2017_samepos,HGMD2017_flank5,avsnp150,gnomad211_genome,gnomad211_exome,exac03,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_afr,1000g2015aug_sas,dbnsfp42a,dbscsnv11,mitotip,GBFreq,gnomADchrM312 -operation f,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --otherinfo --nastring . --vcfinput -out /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.annoVar /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.vcf /biodata/databases/annovar_humandb"
$docker_cmd "perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/gatk4_germline_SNVIndel/merge_result.v1.02.pl -exomiser  /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/exomiser.merge_all.variants.tsv /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.annoVar.hg19_multianno.vcf /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.annoVar.hg19_multianno.txt >/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.annoVar.hg19_merge.txt"
$docker_cmd "perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/gatk4_germline_SNVIndel/filter_clinical_variant.v1.01.pl /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.annoVar.hg19_merge.txt >/mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.annoVar.hg19_merge.filter.txt"
$docker_cmd "cat /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.annoVar.hg19_merge.filter.txt|perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/mongodb/anno2mongodb.v1.01.pl Genetic.Tests A202208020001 "
