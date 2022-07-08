#cat rs_ids.txt |while read L;do grep -i $L /home/liumingming/bio_pipeline/software/AnnoVar/20200317/annovar/humandb/hg19_avsnp150.txt >>rs_ids.txt.xls;done
sed -i 's/rs3064744/rs8175347/g' rs_ids.txt.xls
perl get_database.pl rs_ids.txt.xls raw_data.txt  >database.20200430.xls
