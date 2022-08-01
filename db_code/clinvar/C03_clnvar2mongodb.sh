
bio_dir=/data/backup/home/liumingming/biodata
current_dir=`pwd`


current_dir=`pwd`

docker_cmd="docker run --rm -v /:/mnt -v $bio_dir:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u 1020:1022 liumingming1988/biodocker  bash -c "
cat clinvar_20220416/clinvar_20220416.snpeff.vcf |perl bin/stat_raw_path.v1.pl >clinvar_20220416/clinvar_20220416.snpeff.path_stat.xls
$docker_cmd "cat /mnt/$current_dir/clinvar_20220416/clinvar_20220416.snpeff.path_stat.xls|perl /mnt/$current_dir/bin/clinvar_stat2mongodb.v1.pl Genetic.clinvar_stat"
