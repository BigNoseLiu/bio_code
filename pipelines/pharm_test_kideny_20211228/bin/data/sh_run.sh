#head -1000 variants.tsv | perl code_get_dbsnp_info.pl >variants.1000.rs_info.txt
#head -2000 variants.tsv |tail -1000|perl code_get_dbsnp_info.pl >variants.2000.rs_info.txt
#head -3000 variants.tsv |tail -1000|perl code_get_dbsnp_info.pl >variants.3000.rs_info.txt
#head -4000 variants.tsv |tail -1000|perl code_get_dbsnp_info.pl >variants.4000.rs_info.txt
#head -5000 variants.tsv |tail -1000|perl code_get_dbsnp_info.pl >variants.5000.rs_info.txt
#head -6000 variants.tsv |tail -1000|perl code_get_dbsnp_info.pl >variants.6000.rs_info.txt
#head -7000 variants.tsv |tail -1000|perl code_get_dbsnp_info.pl >variants.7000.rs_info.txt
#cat variants.*rs_info.txt>data_variants.dbsnp_info.txt
#rm -rf variants.*rs_info.txt

#cat data_variants.dbsnp_info.txt |perl code_convert_variant.pl >data_rs.GRCh38.haplotypes.tsv
perl code_get_data.pl data_all_kidney_target_drug.txt temp
#cat temp.var.xls|awk -F "\t" '{print $5"\t"$6"\t"$4"\t"$8"\t"$9}'|grep '^chr'|sort|uniq >data_kidney.target_mut.xls
#grep -v REF temp.var.xls | perl code_get_grch38_fa.pl > data_kidney.target_mut.fa


#cat  temp.var.xls |awk -F "\t" '{print $5"\t"$6"\t"$7"\t"$8"\t"$9}'|grep -v REFERENCE|sort|uniq|grep chr >temp.var.avinput
#perl table_annovar.pl  -buildver hg38 -protocol gnomad211_genome -operation f temp.var.avinput humandb/ --remove --otherinfo --nastring . -out data_mut







wc -l temp*
