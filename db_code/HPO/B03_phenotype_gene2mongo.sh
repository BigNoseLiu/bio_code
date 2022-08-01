current_dir=`pwd`
current_version=$current_dir_temp/data_20220422
docker_cmd="docker run --rm -v /:/mnt -v /data/backup/home/liumingming/biodata:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u 1020:1022 liumingming1988/biodocker  bash -c "
$docker_cmd "perl /mnt/$current_dir/bin/hpo_gene2mongodb.v1.pl Genetic.phenotype2gene /mnt/$current_dir/$current_version/phenotype_to_genes.txt /mnt/$current_dir/$current_version/genes_to_phenotype.txt"
