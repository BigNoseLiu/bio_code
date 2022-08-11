#1. input files
sample_id=$1
read_structure1=$2

#2. system configs
ug_id=1000:1000
docker_name=etal/cnvkit
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"
ref=/mnt/3634AD525B23933E//biodata/database/cnvkit_db/hg19/GRch37.refFlat.txt

$docker_cmd "cnvkit.py target /mnt/$1 --annotate /mnt/$ref --split --output /mnt/$2"
