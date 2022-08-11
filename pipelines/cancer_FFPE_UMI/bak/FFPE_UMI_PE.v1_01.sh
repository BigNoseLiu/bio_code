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
biodata_dir=/mnt/3634AD525B23933E/biodata
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata -u $ug_id $docker_name bash -c"

factera=/biodata/software/factera/factera.pl
factera_exon=/biodata/software/factera/exons.noChr_noMT.bed
factera_exonlevel=/biodata/software/factera/exons.noChr_noMT.exonlevel.bed
ref=/biodata/database/gatk_bundle/human_g1k_v37_decoy.fasta
bed_file=/biodata/pipelines/cancer_FFPE_UMI/bed/cancer120.bed
interval_list=/biodata/pipelines/cancer_FFPE_UMI/bed/cancer120.bed.interval_list
gatk_dbsnp=/biodata/database/gatk_bundle/dbsnp_138.b37.vcf
gatk_Mills_and_1000G_gold_standard_indels=/biodata/database/gatk_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf
mutect_gnomad=/biodata/database/gatk_bundle/af-only-gnomad.raw.sites.b37.vcf.gz
cmd_pon=""

#call consensus reads
mkdir -p $out_dir/$sample_id/
out_prefix="/mnt/$out_dir/$sample_id/$sample_id"
unmap_bam="$out_prefix.00.unmap.bam"
umi_unmap_bam="$out_prefix.01.unmap.withUMI.bam"
umi_unmap_fq="$out_prefix.02.unmap.withUMI.fq"
umi_map_sam="$out_prefix.03.map.withUMI.sam"
merge_bam="$out_prefix.04.merge.bam"
markdup_bam="$out_prefix.21.merge.markdup.bam"
markdup_metric="$out_prefix.21.merge.markdup.metrics"
markdup_umi_metric="$out_prefix.21.merge.umi.metrics"
sorted_bam="$out_prefix.21.merge.markdup.sorted.bam"
sorted_multimetrics="$out_prefix.21.merge.markdup.sorted.multiple_metrics"
sorted_HSmetrics="$out_prefix.21.merge.markdup.sorted.HS_metrics"
recal_table="$out_prefix.21.merge.markdup.sorted.recal_table"
recal_bam="$out_prefix.21.merge.markdup.sorted.recal.bam"
somatic_raw_vcf="$out_prefix.21.merge.markdup.sorted.recal.somatic_raw.vcf"
somatic_f1r2="$out_prefix.21.merge.markdup.sorted.recal.somatic_f1r2.tar.gz"
somatic_bam="$out_prefix.21.merge.markdup.sorted.recal.somatic_out.bam"
vardict_out2="$out_prefix.21.consensus_merge.clip.vardict.out"
vardict_vcf2="$out_prefix.21.consensus_merge.clip.vardict.vcf"

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

ref_2bit=/biodata/database/gatk_bundle/human_g1k_v37_decoy.2bit
$docker_cmd "mkdir -p /mnt/$out_dir/$sample_id/factera;perl $factera -o /mnt/$out_dir/$sample_id/factera $sorted_bam $factera_exon $ref_2bit $bed_file"
exit
#fgbio to mark umi
$docker_cmd "picard -Xmx8G FastqToSam F1=/mnt/$raw_fq1 F2=/mnt/$raw_fq2 O=$unmap_bam SAMPLE_NAME=$sample_id"
$docker_cmd "fgbio  -Xmx8G ExtractUmisFromBam --input=$unmap_bam --output=$umi_unmap_bam --read-structure=$read_structure1 $read_structure2 --molecular-index-tags=ZA ZB --single-tag=RX"
$docker_cmd "picard -Xmx8G SamToFastq I=$umi_unmap_bam F=$umi_unmap_fq INTERLEAVE=true"
$docker_cmd "bwa mem -p -t 5 -M $ref $umi_unmap_fq >$umi_map_sam"
$docker_cmd "picard -Xmx8G MergeBamAlignment UNMAPPED=$umi_unmap_bam ALIGNED=$umi_map_sam O=$merge_bam R=$ref SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"

#preprocess for gatk
$docker_cmd "picard -Xmx8G UmiAwareMarkDuplicatesWithMateCigar INPUT=$merge_bam OUTPUT=$markdup_bam METRICS_FILE=$markdup_metric UMI_METRICS_FILE=$markdup_umi_metric OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500  ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
$docker_cmd "picard -Xmx8G AddOrReplaceReadGroups I=$markdup_bam O=$sorted_bam SORT_ORDER=coordinate RGID=$sample_id  RGLB=lb_$sample_id RGPL=illumina RGPU=pu_$sample_id RGSM=$sample_id"
$docker_cmd "gatk BaseRecalibrator -R $ref -I $sorted_bam -O $recal_table --known-sites $gatk_dbsnp --known-sites $gatk_Mills_and_1000G_gold_standard_indels --use-original-qualities"
$docker_cmd "picard -Xmx8G CollectMultipleMetrics I=$sorted_bam O=$sorted_multimetrics PROGRAM=CollectQualityYieldMetrics R=$ref INTERVALS=$interval_list"
$docker_cmd "picard -Xmx8G CollectHsMetrics I=$sorted_bam O=$sorted_HSmetrics R=$ref BAIT_INTERVALS=$interval_list TARGET_INTERVALS=$interval_list"
#call somatic variant with mutect2
$docker_cmd "gatk ApplyBQSR -R $ref -I $sorted_bam -bqsr $recal_table -O $recal_bam  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30  --add-output-sam-program-record --create-output-bam-md5  --use-original-qualities"
$docker_cmd "gatk Mutect2 -R $ref -I $recal_bam -tumor $sample_id $cmd_pon --germline-resource $mutect_gnomad --af-of-alleles-not-in-resource 0.0000025 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -L $interval_list -O $somatic_raw_vcf --f1r2-tar-gz $somatic_f1r2 -bamout $somatic_bam"

exit

#call somatic variant with vardict
in_bam=`echo $recal_bam|sed 's/^\/mnt//'`
out_vardict=`echo $vardict_out2|sed 's/^\/mnt//'`
vardict-java -f 0.01 -z -c 1 -S 2 -E 3 -g 4 -th 4 -N $sample_id -b $in_bam -G /media/rain/d/$ref /media/rain/d/$bed_file >$out_vardict
$docker_cmd "cat $vardict_out2|teststrandbias.R|var2vcf_valid.pl -N $sample_id -E -f 0.01 >$vardict_vcf2"

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
