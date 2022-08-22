#1. input files
bin_dir=$(dirname $0)
ug_id=$1
sample_id=$2
out_dir=$3
raw_fq1=$4
raw_fq2=$5

#2. system configs
docker_name=liumingming1988/biodocker
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"

biodata_dir=$bin_dir/../../..
db=$biodata_dir/databases/kraken2/HPV_db/

mkdir -p $out_dir
$docker_cmd "kraken2 --db $db --threads 10 --report /mnt/$out_dir/$sample_id.report --output /mnt/$out_dir/$sample_id.output --paired /mnt/$raw_fq1 /mnt/$raw_fq2"
