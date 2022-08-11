#1. input files
sample_id=$1
read_structure1=$2
raw_fq1=$3
read_structure2=$4
raw_fq2=$5
out_dir=$6

#2. system configs
ug_id=1000:1000
docker_name=liumingming1988/biodocker
biodata_dir=/media/rain/d/biodata/
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata -u $ug_id $docker_name bash -c"

ref=/biodata/database/gatk_bundle/human_g1k_v37_decoy.fasta
bed_file=/biodata/pipelines/cancer_FFPE_UMI/bed/cancer120.bed
interval_list=/biodata/pipelines/cancer_FFPE_UMI/bed/cancer120.bed.interval_list

#call consensus reads
mkdir -p $out_dir/$sample_id/
out_prefix="/mnt/$out_dir/$sample_id/$sample_id"
unmap_bam="$out_prefix.00.unmap.bam"
umi_unmap_bam="$out_prefix.01.unmap.withUMI.bam"
umi_unmap_fq="$out_prefix.02.unmap.withUMI.fq"
umi_map_sam="$out_prefix.03.map.withUMI.sam"
merge_bam="$out_prefix.04.merge.bam"
group_bam="$out_prefix.05.group.bam"
duplex_metrics="$out_prefix.05.group.duplex_metrics"
consensus_unmap_bam="$out_prefix.06.consensus_unmap.bam"
consensus_filter_bam="$out_prefix.07.consensus_unmap.filter.bam"
consensus_filter_sort_bam="$out_prefix.07.consensus_unmap.filter.sort.bam"
consensus_filter_fq="$out_prefix.08.consensus_unmap.filter.fq"
consensus_map_sam="$out_prefix.09.consensus_unmap.filter.sam"
consensus_merge_bam="$out_prefix.10.consensus_merge.bam"
consensus_merge_metrics="$out_prefix.10.consensus_merge.multiple_metrics"
consensus_merge_HSmetrics="$out_prefix.10.consensus_merge.HS_metrics"
consensus_clip_bam="$out_prefix.11.consensus_merge.clip.bam"
vardict_out="$out_prefix.12.consensus_merge.clip.vardict.out"
vardict_vcf="$out_prefix.12.consensus_merge.clip.vardict.vcf"

$docker_cmd "picard -Xmx8G CollectHsMetrics I=$consensus_merge_bam O=$consensus_merge_HSmetrics R=$ref BAIT_INTERVALS=$interval_list TARGET_INTERVALS=$interval_list"
exit
$docker_cmd "picard -Xmx8G FastqToSam F1=/mnt/$raw_fq1 F2=/mnt/$raw_fq2 O=$unmap_bam SAMPLE_NAME=$sample_id"
$docker_cmd "fgbio  -Xmx8G ExtractUmisFromBam --input=$unmap_bam --output=$umi_unmap_bam --read-structure=$read_structure1 $read_structure2 --molecular-index-tags=ZA ZB --single-tag=RX"
$docker_cmd "picard -Xmx8G SamToFastq I=$umi_unmap_bam F=$umi_unmap_fq INTERLEAVE=true"
$docker_cmd "bwa mem -p -t 5 -M $ref $umi_unmap_fq >$umi_map_sam"
$docker_cmd "picard -Xmx8G MergeBamAlignment UNMAPPED=$umi_unmap_bam ALIGNED=$umi_map_sam O=$merge_bam R=$ref SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
$docker_cmd "fgbio -Xmx8G GroupReadsByUmi --input=$merge_bam --output=$group_bam --strategy=paired --edits=1 --min-map-q=20 --raw-tag=RX"
$docker_cmd "fgbio -Xmx8G CollectDuplexSeqMetrics --input=$group_bam --output=$duplex_metrics --intervals=$interval_list --duplex-umi-counts=true"
$docker_cmd "fgbio -Xmx8G CallDuplexConsensusReads --input=$group_bam --output=$consensus_unmap_bam --error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=20 --min-reads=1"
$docker_cmd "fgbio -Xmx8G FilterConsensusReads --input=$consensus_unmap_bam --output=$consensus_filter_bam --ref=$ref --min-reads=2 1 1 --max-read-error-rate=0.05 --max-base-error-rate=0.1 --min-base-quality=50 --max-no-call-fraction=0.05"
$docker_cmd "picard -Xmx8G SamToFastq I=$consensus_filter_bam F=$consensus_filter_fq INTERLEAVE=true"
$docker_cmd "bwa mem -p -t 5 -M $ref $consensus_filter_fq >$consensus_map_sam"

$docker_cmd "picard -Xmx8G SortSam I=$consensus_filter_bam O=$consensus_filter_sort_bam SO=queryname"
$docker_cmd "picard -Xmx8G MergeBamAlignment UNMAPPED=$consensus_filter_sort_bam ALIGNED=$consensus_map_sam O=$consensus_merge_bam R=$ref SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true ATTRIBUTES_TO_RETAIN=XS"
$docker_cmd "picard -Xmx8G CollectMultipleMetrics I=$consensus_merge_bam O=$consensus_merge_metrics R=$ref INTERVALS=$interval_list"
$docker_cmd "picard -Xmx8G CollectHsMetrics I=$consensus_merge_bam O=$consensus_merge_HSmetrics R=$ref BAIT_INTERVALS=$interval_list TARGET_INTERVALS=$interval_list"
##$docker_cmd "fgbio -Xmx8G ClipBAM --input=$consensus_merge_bam --output=$consensus_clip_bam --ref=$ref --soft-clip=false --clip-overlapping-reads=true"

##$docker_cmd "vardict -f 0.01 -z -c 1 -S 2 -E 3 -g 4 -th 4 -N $sample_id -b $consensus_merge_bam -G $ref $bed_file|teststrandbias.R|var2vcf_valid.pl -N $sample_id -E -f 0.01 >$vardict_vcf"
in_bam=`echo $consensus_merge_bam|sed 's/^\/mnt//'`
out_vardict=`echo $vardict_out|sed 's/^\/mnt//'`
vardict-java -f 0.01 -z -c 1 -S 2 -E 3 -g 4 -th 4 -N $sample_id -b $in_bam -G /media/rain/d/$ref /media/rain/d/$bed_file >$out_vardict
$docker_cmd "cat $vardict_out|teststrandbias.R|var2vcf_valid.pl -N $sample_id -E -f 0.01 >$vardict_vcf"
exit



docker run --rm -v $cw:/mnt -u 1001:1001 liumm_docker bash -c "picard MarkDuplicates I=/mnt/work/kefu_20211222-ChangShaFuYinR/align.sam O=/mnt/work/kefu_20211222-ChangShaFuYinR/align.sam.markDup.bam M=/mnt/work/kefu_20211222-ChangShaFuYinR/align.sam.markDup.metric VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname"
docker run --rm -v $cw:/mnt -u 1001:1001 liumm_docker bash -c "picard AddOrReplaceReadGroups I=/mnt/work/kefu_20211222-ChangShaFuYinR/align.sam.markDup.bam O=/mnt/work/kefu_20211222-ChangShaFuYinR/align.sam.markDup.sort.bam SORT_ORDER=coordinate RGID=sample_id  RGLB=lb_sample_id RGPL=illumina RGPU=pu_sample_id RGSM=sample_id TMP_DIR=$cw"
docker run --rm -v $cw:/mnt -u 1001:1001 liumm_docker bash -c "samtools index /mnt/work/kefu_20211222-ChangShaFuYinR/align.sam.markDup.sort.bam"
