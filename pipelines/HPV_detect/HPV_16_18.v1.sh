#1. input files
sample_id=$1
raw_fq1=$2
raw_fq2=$3
out_dir=$4

#2. system configs
ug_id=0:0
docker_name=liumingming1988/biodocker
biodata_dir=/home/rain/liumm_data/biodata
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata -u $ug_id $docker_name bash -c"

db=/biodata/databases/kraken2/HPV_db/
mkdir -p $out_dir
$docker_cmd "kraken2 --db $db --threads 10 --report /mnt/$out_dir/$sample_id.report --output /mnt/$out_dir/$sample_id.output --paired /mnt/$raw_fq1 /mnt/$raw_fq2"
