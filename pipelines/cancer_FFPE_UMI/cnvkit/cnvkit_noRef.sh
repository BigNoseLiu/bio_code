#1. input files
bin_dir=$(dirname $0)
ug_id=$1
sample_id=$2
out_dir=$3
tumor_dir=$4
normal_dir=$5

cwd=`pwd`
#2. system configs
docker_name=etal/cnvkit
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"
ref_flat=$bin_dir/hg19.refFlat.no_chr.txt
bed=$bin_dir/cancer120.cnvkit.bed
ref_fa=$bin_dir/../../../../databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta
ref_access=$bin_dir/access-5k-mappable.grch37.bed;

$docker_cmd "mkdir -p /mnt/$out_dir/$sample_id"
$docker_cmd "cnvkit.py batch /mnt/$tumor_dir/*.bam --normal /mnt/$normal_dir/*.bam --targets /mnt/$bed --annotate /mnt/$ref_flat --fasta /mnt/$ref_fa --access /mnt/$ref_access  --output-reference /mnt/$out_dir/reference.cnn --output-dir /mnt/$out_dir/$sample_id --diagram --scatter -p 10"
