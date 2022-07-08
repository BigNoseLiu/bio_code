#处理bed文件,note: interval必须为interval_list结尾
#sed 's/chr//' Exome-Agilent_V6.bed >Exome-Agilent_V6.no_chr.bed
#java -jar ../../../software/Picard/20180930/picard.jar BedToIntervalList I=Exome-Agilent_V6.no_chr.bed O=Exome-Agilent_V6.no_chr.bed.interval_list SD=../../../database/genomes/b37_decoy/human_g1k_v37_decoy.dict

sh fix4gatk.sh yubing_cancer/NGS
exit
sh fix4gatk.sh cancer600gene/cancer600gene
sh fix4gatk.sh zhongkeV6/Aglient_V6
sh fix4gatk.sh Exome-Agilent_V6/Exome-Agilent_V6
sh fix4gatk.sh Exome-Agilent_SureSelect_V7/S31285117_hs_hg19/S31285117_Regions
#bed_file=lianchuan/cancer_gene600/lianchuan_600gene
#cat $bed_file\.bed|perl fix_bed_format_for_gatk.pl >$bed_file.fix4gatk.bed
#java -jar ../../../software/Picard/20180930/picard.jar BedToIntervalList I=$bed_file.fix4gatk.bed O=$bed_file.fix4gatk.bed.interval_list SD=../../../database/genomes/b37_decoy/human_g1k_v37_decoy.dict
#exit
#sed 's/chr//' $bed_file\.bed >$bed_file.no_chr.bed
#java -jar ../../../software/Picard/20180930/picard.jar BedToIntervalList I=$bed_file.no_chr.bed O=$bed_file.no_chr.bed.interval_list SD=../../../database/genomes/b37_decoy/human_g1k_v37_decoy.dict
