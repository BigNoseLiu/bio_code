current_dir=`pwd`
current_version=$current_dir_temp/data_20220422
docker_cmd="docker run --rm -v /:/mnt -v /data/backup/home/liumingming/biodata:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u 1020:1022 liumingming1988/biodocker  bash -c "
#$docker_cmd "perl /mnt/$current_dir/../bin/delete_all_mongodb.pl Genetic.phenotypes"
$docker_cmd "cat /mnt/$current_dir/$current_version/hp.obo|perl /mnt/$current_dir/bin/hpo2mongodb.v2.pl Genetic.phenotypes /mnt/$current_dir/chpo.201802.txt /mnt/$current_dir/$current_version/phenotype_to_genes.txt /mnt/$current_dir/$current_version/genes_to_phenotype.txt"
