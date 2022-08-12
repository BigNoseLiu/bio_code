#1. input files
cwd=`pwd`
#2. system configs
ug_id=1000:1000
docker_name=etal/cnvkit
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"
ref_flat=/mnt/3634AD525B23933E//biodata/database/cnvkit_db/hg19/GRch37.refFlat.txt
bed=/mnt/3634AD525B23933E/biodata/pipelines/cancer_FFPE_UMI/bed/cancer120.cnvkit.bed
ref_fa=/mnt/3634AD525B23933E/biodata/database/gatk_bundle/human_g1k_v37_decoy.fasta
ref_access=/mnt/3634AD525B23933E/biodata/database/gatk_bundle/access-5kb.human_g1k_v37_decoy.bed

$docker_cmd "cnvkit.py batch /mnt/$cwd/tumor_samples/*.bam --normal /mnt/$cwd/normal_samples/*.bam --targets /mnt/$bed --annotate /mnt/$ref_flat --fasta /mnt/$ref_fa --access /mnt/$ref_access  --output-reference /mnt/$cwd/cancer120_reference.cnn --output-dir /mnt/$cwd/tumor_samples --diagram --scatter -p 3"
