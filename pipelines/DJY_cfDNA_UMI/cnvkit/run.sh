ug_id=1000:1000
docker_name=etal/cnvkit
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"

cwd=`pwd`
bed=$cwd/../DJY_bed/DJY_predict.bed
out_bed=$cwd/DJY_predict.cnvkit.bed
ref_flat=$cwd/hg19.refFlat.no_chr.txt


$docker_cmd "cnvkit.py target /mnt/$bed --annotate /mnt/$ref_flat --split --output /mnt/$out_bed"
