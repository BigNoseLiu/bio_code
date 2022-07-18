bio_dir=/data/backup/home/liumingming/biodata
version=20220416

current_dir=$bio_dir/databases/clinvar/clinvar_$version
docker_cmd="docker run --rm -v /:/mnt -v /data/backup/home/liumingming/biodata:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u 1020:1022 liumingming1988/biodocker  bash -c "
$docker_cmd "java -jar /biodata/software/snpeff/snpEff5_0/snpEff.jar -v GRCh37.p13.RefSeq -s /mnt/$current_dir/clinvar_20220416.snpeff.stat.html /mnt/$current_dir/clinvar_20220416.vcf.gz >/mnt/$current_dir/clinvar_20220416.snpeff.vcf"
