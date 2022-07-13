for db in 'dbnsfp42a' 'avsnp150' 'gnomad211_exome' 'gnomad30_genome' 'dbscsnv11' 'revel' '1000g2015aug'
do
	perl bin/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar $db ../../../databases/annovar_humandb/
done
