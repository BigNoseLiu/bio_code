#1. input files
bin_dir=$(dirname $0)
ug_id=$1
sample_id=$2
out_dir=$3
raw_fq1=$4
raw_fq2=$5

#2. system configs
docker_name=rachelbj/hla-hd
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"

#fgbio to mark umi
$docker_cmd "mkdir -p /mnt/$out_dir"
$docker_cmd "/app/hlahd.1.4.0/bin/hlahd.sh -t 6 -m 100 -c 0.95 -f /app/hlahd.1.4.0/freq_data /mnt/$raw_fq1 /mnt/$raw_fq2 /app/hlahd.1.4.0/HLA_gene.split.3.32.0.txt /app/hlahd.1.4.0/dictionary $sample_id /mnt/$out_dir"
