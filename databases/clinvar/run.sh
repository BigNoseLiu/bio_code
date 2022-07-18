current_dir=`pwd`
docker_cmd="docker run --rm -v /:/mnt -v /data/backup/home/liumingming/biodata:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u 1020:1022 liumingming1988/biodocker  bash -c "
$docker_cmd "java -jar /biodata/software/snpeff/snpEff5_0/snpEff.jar -v GRCh37.p13.RefSeq -s /mnt/$current_dir/clinvar_20220416.snpeff.stat.html /mnt/$current_dir/clinvar_20220416.vcf.gz >/mnt/$current_dir/clinvar_20220416.snpeff.vcf"
cat clinvar_20220416.snpeff.vcf |perl ../bin/stat_raw_path.v1.pl >clinvar_20220416.snpeff.path_stat.xls
$docker_cmd "perl /mnt/$current_dir/../bin/delete_all_mongodb.pl Genetic.clinvar_stat"
$docker_cmd "cat /mnt/$current_dir/clinvar_20220416.snpeff.path_stat.xls|perl /mnt/$current_dir/../bin/clinvar_stat2mongodb.v1.pl Genetic.clinvar_stat"
exit
#perl /data/backup/home/liumingming/biodata/biopipe/human_germline/anno_vcf/bin_202204/vcf_tools/fix_vcf4_2.v1.pl -in clinvar_20220416.vcf -out clinvar_20220416.fix.vcf -ref /data/backup/home/liumingming/biodata/databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta
#$docker_cmd "bgzip -f /mnt/$current_dir/clinvar_20220416.fix.vcf"
$docker_cmd "tabix -p vcf /mnt/$current_dir/clinvar_20220416.fix.vcf.gz"
$docker_cmd "bcftools norm -f /biodata/databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta -o /mnt/$current_dir/clinvar_20220416.fix.bcftools_leftAlign.vcf /mnt/$current_dir/clinvar_20220416.fix.vcf.gz"
$docker_cmd "java -jar /biodata/software/snpeff/snpEff5_0/snpEff.jar -v GRCh37.p13.RefSeq -s /mnt/$current_dir/clinvar_20220416.snpeff.stat.html /mnt/$current_dir/clinvar_20220416.fix.bcftools_leftAlign.vcf >/mnt/$current_dir/clinvar_20220416.fix.bcftools_leftAlign.snpeff.vcf"



#prepare for annovar
perl ../../../software/annovar/annovar/convert2annovar.pl -format vcf4  clinvar_20220416.vcf.gz -outfile clinvar_20220416.conver2annovar -includeinfo
echo -e "#Chr\tStart\tEnd\tRef\tAlt\tclinvar_anno" >hg19_clinvar_20220416.txt
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\tPOS="$1":"$2":"$3":"$4":"$5";CLNVARID="$8";"$13}' clinvar_20220416.conver2annovar|perl ../bin/fix_comma.pl >>hg19_clinvar_20220416.txt
perl ../../../software/annovar/buid_index.pl -in hg19_clinvar_20220416.txt -bin 1000

sed '1d' hg19_clinvar_20220416.txt |awk -F "\t" '{print "1\t"$1"\t"$2-5"\t"$3+5"\t"$6}' >hg19_clinvar_region_20220416.txt

