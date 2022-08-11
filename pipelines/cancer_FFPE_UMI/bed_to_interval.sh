#1. input files
sample_id=$1
read_structure1=$2

#2. system configs
ug_id=1000:1000
docker_name=liumingming1988/biodocker
biodata_dir=/media/rain/d/biodata/
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata -u $ug_id $docker_name bash -c"
ref=/biodata/database/gatk_bundle/human_g1k_v37_decoy.fasta

$docker_cmd "picard -Xmx8G BedToIntervalList I=/mnt/$1 O=/mnt/$2 SD=$ref"
