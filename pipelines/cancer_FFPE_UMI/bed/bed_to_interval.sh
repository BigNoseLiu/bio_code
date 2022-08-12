#1. input files
bin_dir=$(dirname $0)

#2. system configs
ug_id=1000:1000
docker_name=liumingming1988/biodocker
docker_cmd="docker run --rm -v /:/mnt -v $bin_dir/../../../../:/biodata -u $ug_id $docker_name bash -c"
ref=/biodata/databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta

$docker_cmd "picard -Xmx8G BedToIntervalList I=/mnt/$1 O=/mnt/$2 SD=$ref"
