#2. system configs
ug_id=1000:1000
docker_name=liumingming1988/biodocker
biodata_dir=/mnt/3634AD525B23933E/biodata/
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata -u $ug_id $docker_name bash -c"
cwd=`pwd`
gunzip human_g1k_v37_decoy.fasta.gz
$docker_cmd " faToTwoBit /mnt/$cwd/human_g1k_v37_decoy.fasta /mnt/$cwd/human_g1k_v37_decoy.2bit"
$docker_cmd "bwa index /mnt/$cwd/human_g1k_v37_decoy.fasta"
