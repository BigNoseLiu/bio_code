bio_dir=/data/backup/home/liumingming/biodata
version=20220416

annovar_db_dir=$bio_dir/databases/annovar_humandb/
clinvar_dir=$bio_dir/databases/clinvar/clinvar_$version
sed '1d' $annovar_db_dir/hg19_clinvar_$version.txt |awk -F "\t" '{print "1\t"$1"\t"$2"\t"$3"\t"$6}' >$annovar_db_dir/hg19_clinvar_samepos_$version.txt
exit
#prepare for annovar
#perl ../../software/annovar/bin/convert2annovar.pl -format vcf4  $clinvar_dir/clinvar_$version.vcf.gz -outfile $clinvar_dir/clinvar_$version.conver2annovar -includeinfo
echo -e "#Chr\tStart\tEnd\tRef\tAlt\tclinvar_anno" >$annovar_db_dir/hg19_clinvar_$version.txt
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\tPOS="$1":"$2":"$3":"$4":"$5";CLNVARID="$8";"$13}' $clinvar_dir/clinvar_20220416.conver2annovar|perl bin/fix_comma.pl >>$annovar_db_dir/hg19_clinvar_$version.txt
perl ../../software/annovar/bin/buid_index.pl -in $annovar_db_dir/hg19_clinvar_$version.txt -bin 1000

sed '1d' $annovar_db_dir/hg19_clinvar_$version.txt |awk -F "\t" '{print "1\t"$1"\t"$2-5"\t"$3+5"\t"$6}' >$annovar_db_dir/hg19_clinvar_flank5_$version.txt
sed '1d' $annovar_db_dir/hg19_clinvar_$version.txt |awk -F "\t" '{print "1\t"$1"\t"$2"\t"$3"\t"$6}' >$annovar_db_dir/hg19_clinvar_samepos_$version.txt

