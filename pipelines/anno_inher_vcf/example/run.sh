ug_id=1020:1022
docker_name=liumingming1988/biodocker
biodata_dir=/data/backup/home/liumingming/biodata/
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u $ug_id $docker_name bash -c "
current_dir=`pwd`
in_vcf=../../human_germline/anno_vcf/example/out_data/germline_anno/test/test.final.leftalign.annoGene.annoVar.hg19_multianno.txt
#head -1 $in_vcf >test.annovar.xls
#cat $in_vcf|awk -F "\t" '$13<=0.05' >>test.annovar.xls
$docker_cmd "perl /mnt/$current_dir/../bin_v2/delete_all_mongodb.pl Genetic.Tests"
#$docker_cmd "cat /mnt/$current_dir/test.annovar.xls|perl /mnt/$current_dir/../bin_v2/anno2mongodb.v1.pl Genetic.Tests A10010101 "
