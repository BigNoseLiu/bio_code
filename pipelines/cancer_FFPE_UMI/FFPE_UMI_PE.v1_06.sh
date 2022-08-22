#1. input files
bin_dir=$(dirname $0)
#ug_id=$1
ug_id=`id|awk -F "[= (]" '{print $2":"$5}'`
sample_id=$2
out_dir=$3
read_structure1=$4
raw_fq1=$5
read_structure2=$6
raw_fq2=$7

#2. system configs
docker_name=liumingming1988/biodocker
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"

biodata_dir=$bin_dir/../../..
annovar=$biodata_dir/git_code/software/annovar/bin/table_annovar.pl
humandb=$biodata_dir/databases/annovar_humandb
factera=$biodata_dir/git_code/software/factera/factera.pl
factera_exon=$biodata_dir/git_code/software/factera/exons.no_chr.bed
ref=$biodata_dir/databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta
ref_2bit=$biodata_dir/databases/gatk_bundle/b37/human_g1k_v37_decoy.2bit
bed_file=$biodata_dir/git_code/pipelines/cancer_FFPE_UMI/bed/cancer120.bed
interval_list=$biodata_dir/git_code/pipelines/cancer_FFPE_UMI/bed/cancer120.bed.interval_list
gatk_dbsnp=$biodata_dir/databases/gatk_bundle/b37/dbsnp_138.b37.vcf
gatk_Mills_and_1000G_gold_standard_indels=$biodata_dir/databases/gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
mutect_gnomad=$biodata_dir/databases/gatk_bundle/b37/af-only-gnomad.raw.sites.b37.vcf.gz
gatk_GetPileupSummaries=$biodata_dir/databases/gatk_bundle/b37/small_exac_common_3_b37.vcf.gz
cmd_pon=""

#call consensus reads
mkdir -p $out_dir/$sample_id/ $out_dir/temp_dir
temp_dir="/mnt/$out_dir/temp_dir"
out_prefix="/mnt/$out_dir/$sample_id/$sample_id"
unmap_bam="$out_prefix.00.unmap.bam"
umi_unmap_bam="$out_prefix.01.unmap.withUMI.bam"
umi_unmap_fq="$out_prefix.02.unmap.withUMI.fq"
umi_map_sam="$out_prefix.03.map.withUMI.sam"
merge_bam="$out_prefix.04.merge.bam"
markdup_bam="$out_prefix.05.markdup.bam"
markdup_metric="$out_prefix.05.markdup.metrics"
markdup_umi_metric="$out_prefix.05.markdup_umi.metrics"
sorted_bam="$out_prefix.06.sorted.bam"
sorted_multimetrics="$out_prefix.07.multiple_metrics"
sorted_HSmetrics="$out_prefix.07.HS_metrics"
recal_table="$out_prefix.08.recal_table"
recal_bam="$out_prefix.08.recal.bam"
recal_multimetrics="$out_prefix.08.multiple_metrics"
recal_HSmetrics="$out_prefix.08.HS_metrics"
somatic_raw_vcf="$out_prefix.09.somatic_raw.vcf"
somatic_f1r2="$out_prefix.09.somatic_f1r2.tar.gz"
table_getpileupsummaries="$out_prefix.09.getpileupsummaries.table"
table_calculatecontamination="$out_prefix.09.calculatecontamination.table"
table_segmentation="$out_prefix.09.segments.table"
somatic_rom="$out_prefix.09.read-orientation-model.tar.gz"
somatic_bam="$out_prefix.09.somatic_out.bam"
somatic_raw_annovar="$out_prefix.09.somatic_raw.annovar"
somatic_oncefiltered_vcf="$out_prefix.10.somatic_raw.filter.vcf"
somatic_oncefiltered_annovar="$out_prefix.10.somatic_raw.filter.annovar"

#fgbio to mark umi
$docker_cmd "picard -Xmx8G FastqToSam F1=/mnt/$raw_fq1 F2=/mnt/$raw_fq2 O=$unmap_bam SAMPLE_NAME=$sample_id"
$docker_cmd "fgbio  -Xmx8G ExtractUmisFromBam --input=$unmap_bam --output=$umi_unmap_bam --read-structure=$read_structure1 $read_structure2 --molecular-index-tags=ZA ZB --single-tag=RX"
$docker_cmd "picard -Xmx8G SamToFastq I=$umi_unmap_bam F=$umi_unmap_fq INTERLEAVE=true"
$docker_cmd "bwa mem -p -t 5 -M $ref $umi_unmap_fq >$umi_map_sam"
$docker_cmd "picard -Xmx8G MergeBamAlignment UNMAPPED=$umi_unmap_bam ALIGNED=$umi_map_sam O=$merge_bam R=$ref SO=coordinate ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true  MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant UNMAPPED_READ_STRATEGY=COPY_TO_TAG ALIGNER_PROPER_PAIR_FLAGS=true UNMAP_CONTAMINANT_READS=true TMP_DIR=$temp_dir"

#preprocess for gatk
$docker_cmd "picard -Xmx8G UmiAwareMarkDuplicatesWithMateCigar INPUT=$merge_bam OUTPUT=$markdup_bam METRICS_FILE=$markdup_metric UMI_METRICS_FILE=$markdup_umi_metric OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500  ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true  TMP_DIR=$temp_dir"
$docker_cmd "picard -Xmx8G AddOrReplaceReadGroups I=$markdup_bam O=$sorted_bam SORT_ORDER=coordinate RGID=$sample_id  RGLB=lb_$sample_id RGPL=illumina RGPU=pu_$sample_id RGSM=$sample_id  TMP_DIR=$temp_dir"
$docker_cmd "gatk BaseRecalibrator -R $ref -I $sorted_bam -O $recal_table --known-sites $gatk_dbsnp --known-sites $gatk_Mills_and_1000G_gold_standard_indels --use-original-qualities"
$docker_cmd "gatk ApplyBQSR -R $ref -I $sorted_bam -bqsr $recal_table -O $recal_bam  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30  --add-output-sam-program-record --create-output-bam-md5  --use-original-qualities"
#QC metrics
$docker_cmd "picard -Xmx8G CollectMultipleMetrics I=$sorted_bam O=$sorted_multimetrics PROGRAM=CollectQualityYieldMetrics R=$ref INTERVALS=$interval_list  TMP_DIR=$temp_dir"
$docker_cmd "picard -Xmx8G CollectHsMetrics I=$sorted_bam O=$sorted_HSmetrics R=$ref BAIT_INTERVALS=$interval_list TARGET_INTERVALS=$interval_list  TMP_DIR=$temp_dir"
$docker_cmd "picard -Xmx8G CollectMultipleMetrics I=$recal_bam O=$recal_multimetrics PROGRAM=CollectQualityYieldMetrics R=$ref INTERVALS=$interval_list  TMP_DIR=$temp_dir"
$docker_cmd "picard -Xmx8G CollectHsMetrics I=$recal_bam O=$recal_HSmetrics R=$ref BAIT_INTERVALS=$interval_list TARGET_INTERVALS=$interval_list  TMP_DIR=$temp_dir"
#call somatic variant with mutect2
$docker_cmd "gatk Mutect2 -R $ref -I $recal_bam -tumor $sample_id $cmd_pon --germline-resource $mutect_gnomad --af-of-alleles-not-in-resource 0.0000025 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -L $interval_list -O $somatic_raw_vcf --f1r2-tar-gz $somatic_f1r2 -bamout $somatic_bam"
#filter mutect2 raw variant
$docker_cmd "gatk LearnReadOrientationModel -I $somatic_f1r2 -O $somatic_rom"
$docker_cmd "gatk GetPileupSummaries -R $ref -V $gatk_GetPileupSummaries --interval-set-rule INTERSECTION -L $interval_list -L $gatk_GetPileupSummaries -I $recal_bam -O $table_getpileupsummaries"
$docker_cmd "gatk CalculateContamination -I $table_getpileupsummaries -O $table_calculatecontamination --tumor-segmentation $table_segmentation"
$docker_cmd "gatk FilterMutectCalls -R $ref -V $somatic_raw_vcf --contamination-table $table_calculatecontamination -O $somatic_oncefiltered_vcf --tumor-segmentation $table_segmentation --ob-priors $somatic_rom"
#annovar variant with annovar
$docker_cmd "perl $annovar -buildver hg19 -protocol refGeneWithVer -operation gx $somatic_oncefiltered_vcf $humandb  --remove --otherinfo --nastring . --vcfinput -out $somatic_oncefiltered_annovar"
