#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(rmtree);
use Cwd qw(abs_path);
use Data::Dumper;
use lib "$Bin/../../../models/perl/lib";
use biopipe;

my ( $help, $infile, $out_dir, $GATK_PON, $optional_cmd, $debug );
GetOptions(
	"help"=>\$help,
	"debug"=>\$debug,
	"opt"=>\$optional_cmd,
	"in=s"=>\$infile,
	"out=s"=>\$out_dir,
	"pon=s"=>\$GATK_PON,
);

#打印帮助文件
if ( defined($help) || !defined($infile) ) {
	&usageWithColor();
	exit 0;
}

#获取程序所在路径
my $bin_dir = "";
my $bio_dir = "";
if( $Bin =~ /^(.*)\/(biodata\/.*)$/ ){
	$bio_dir = "$1/biodata";
	$bin_dir = "/$2";
}
else{
	print STDERR "Err: you're not allowed to change the name of biodata!!!\n";
	exit(0);
}
#设置docker运行命令
if( !defined($out_dir) || $out_dir !~ /\S/ ){
	my $current_dir = `pwd`;
	chomp $current_dir;
	$out_dir = "$current_dir/data_out";
	`mkdir -p $out_dir`;
}
elsif( !(-d $out_dir) ){
	print STDERR "Err: $out_dir not exist or not a directory!!!\n";
	exit(0);
}
#初始化docker相关变量
&init_biopipe( $out_dir, $bio_dir, "li/biojupyter" );


#软件变量配置
my $bin_mkdir		=	"$docker_cmd mkdir -p";
my $bin_fastp		=	"$docker_cmd fastp";
my $bin_trimmomatic	=	"$docker_cmd trimmomatic";
my $bin_polysolver="/usr/bin/docker run -P --name \$NAME -v \$DIR1:\$DIR2 sachet/polysolver:v4 bash /home/polysolver/scripts/shell_call_hla_type \$DIR2/\$BAM Unknown 1 hg19 STDFQ 0 \$DIR2/\$OUT_DIR; /usr/bin/docker rm \$NAME";
my $bin_statHLA="perl $bin_dir/HLA/stat_HLA.v1.01.pl";
#my $bin_polysolver="/usr/bin/docker run -P --name \$NAME -v \$DIR1:\$DIR2 sachet/polysolver:v4 bash /home/polysolver/scripts/shell_call_hla_type \$DIR2/\$BAM Unknown 1 hg19 STDFQ 0 \$DIR2; /usr/bin/docker rm \$NAME";
my $bin_msisensor	=	"$docker_cmd msisensor-pro";
my $bin_factera		=	"$docker_cmd perl /biodata/software/FACTERA/v1.4.4/factera.pl";
my $factera_bed		=	"/biodata/software/FACTERA/v1.4.4/exons.noChr_noMT.bed";
my $factera_bed2	=	"/biodata/software/FACTERA/v1.4.4/exons.noChr_noMT.exonlevel.bed";
my $bin_facterastat	=	"$docker_cmd perl $bin_dir/factera/stat_factera.pl";
my $bin_filter_umi_cfbestV2=	"$docker_cmd perl $bin_dir/umi/trimfq_cfbest_umi_v2.pl";
my $bin_fgbio		=	"$docker_cmd fgbio";
my $bin_varscan		=	"$docker_cmd varscan";
my $bin_snpeff		=	"$docker_cmd java -jar /biodata/software/snpeff/20181017/snpEff/snpEff.jar -v GRCh37.p13.RefSeq";
my $bin_select_snpeff	=	"perl $bin_dir/snpeff/select_standard_hgvs_from_snpeff_Eff.v1.01.pl";
my $bin_anno_gene	=	"perl $bin_dir/var_inter/anno_gene.v2.01.pl";
my $bin_fix_snpeff	=	"perl $bin_dir/snpeff/fix_snpeff_Eff.v1.pl";
#my $bin_fix_snpeff	=	"$docker_cmd perl $bin_dir/snpeff/fix_snpeff_Eff.v1.pl|perl $bin_dir/snpeff/select_standard_hgvs_from_snpeff_Eff.v1.01.pl|perl  $bin_dir/var_inter/anno_gene.v2.01.pl`";
#"cat $snpeff_anno_vcf|$bin_fix_snpeff|$bin_select_snpeff|$bin_anno_gene -clinvar_all $database_clinvar_stat_all -clinvar_mul $database_clinvar_stat_mul -hpo_g2p $database_g2p_hpo  -hpo_p2g $database_p2g_hpo >$snpeff_anno_fix_vcf",
my $bin_bedtools	=	"$docker_cmd bedtools";
my $bin_qualimap="$bin_dir/../../software/QualiMap/qualimap_v2.2.1/qualimap";
my $bin_filterMutect	=	"$docker_cmd perl $bin_dir/gatk4_somatic_SNVIndel/filter_mutect2_pair.leftSplit.pl";
my $bin_fixVcf4		=	"$docker_cmd perl $bin_dir/gatk4_germline_SNVIndel/fix_gatk4_vcf4_2.v2.pl";
my $bin_tag_nearbyOverlap_variant="$docker_cmd perl $bin_dir/gatk4_germline_SNVIndel/tag_nearbyOverlap_variant.v1.01.pl";
my $bin_statsex		=	"$docker_cmd perl $bin_dir/sex_stat/stat_sex.v1.01.pl";
my $bin_statQC		=	"$docker_cmd perl $bin_dir/QC/stat_qc.v1.01.pl";
my $bin_statPicardQC	=	"$docker_cmd perl $bin_dir/QC/stat_picardQC.v1.01.pl";
my $bin_statmsi		=	"$docker_cmd perl $bin_dir/msi/stat_msi.v1.01.pl";
my $bin_smallVariant	=	"$docker_cmd perl $bin_dir/gatk4_germline_SNVIndel/stat_smallVariant.v1.01.pl";
my $bin_bwa		=	"$docker_cmd bwa";
my $bin_picard		=	"$docker_cmd picard";
my $bin_samtools	=	"$docker_cmd samtools";
my $bin_bcftools	=	"$docker_cmd bcftools";
#exomiser
my $bin_exomiser	=	"$docker_cmd java -Xms2g -Xmx4g -jar /biodata/software/exomiser/software/exomiser-cli-12.1.0/exomiser-cli-12.1.0.jar";
my $config_exomiser	=	"/biodata/software/exomiser/software/exomiser-cli-12.1.0/application.properties";
my $bin_create_exomiser_exome=	"$docker_cmd cat $bin_dir/exomiser/exome.yml|perl $bin_dir/exomiser/create_yml.pl";
#annovar
my $bin_annovar		=	"/biodata/software/annovar/20200317/annovar/";
my $bin_annovar_withPara=	"$docker_cmd perl $bin_annovar/table_annovar.pl -buildver hg19 -protocol refGeneWithVer,cosmic70,cosmic81,civic_20200317,civicRegion4SmallMut_20200317,clinvar_20190305,HGMD2017,gwasCatalog,avsnp150,gnomad211_genome,gnomad211_exome,esp6500siv2_all,dbnsfp35a,nci60,1000g2014oct_all,1000g2014oct_eas,1000g2014oct_eur,1000g2014oct_amr,1000g2014oct_afr,1000g2014oct_sas,ljb26_all -operation gx,f,f,f,r,f,f,r,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --otherinfo --nastring . --vcfinput";
my $bin_table_annovar = 	"$docker_cmd perl $bin_annovar/table_annovar.pl";
my $annovar_humandb	=	"$bin_annovar/humandb";
#HPO:gene to disease
my $database_g2p_hpo	=	"/biodata/database/HPO/20200326/genes_to_phenotype.txt";
my $database_p2g_hpo	=	"/biodata/database/HPO/20200326/phenotype_to_genes.txt";
#ClinVar, note: must be matched with $bin_snpeff(eg. -v GRCh37.p13.RefSeq)
my $database_clinvar_stat_all=	"/biodata/database/clinvar/20200326/clinvar_20200316.snpeff.fix.vcf.filter_no.stat.xls";
my $database_clinvar_stat_mul=	"/biodata/database/clinvar/20200326/clinvar_20200316.snpeff.fix.vcf.filter_multiple.stat.xls";
#ref genome
my $ref_fa		=	"/biodata/database/genomes/b37_decoy/human_g1k_v37_decoy.fasta";
my $ref_fa_dict		=	"/biodata/database/genomes/b37_decoy/human_g1k_v37_decoy.dict";
my $ref_fa_2bit		=	"/biodata/database/genomes/b37_decoy/human_g1k_v37_decoy.2bit";
my $ref_flat		=	"/biodata/software/CNVkit/database/hg19/GRch37.refFlat.txt";
my $ref_msi_microsatellite=	"/biodata/database/genomes/b37_decoy/human_g1k_v37_decoy.microsatellites.list";
#varDict
my $bin_vardict="$bin_dir/../../software/VarDict/VarDictJava-master/VarDict/VarDict-1.8.0/bin/VarDict";
my $bin_vardict_teststrandbias="$bin_dir/../../software/VarDict/VarDict-master/teststrandbias.R";
my $bin_vardict_var2vcf_valid="perl $bin_dir/../../software/VarDict/VarDict-master/var2vcf_valid.pl";
#gatk
my $bin_gatk		=	"$docker_cmd gatk";
my $gatk_dbsnp		=	"/biodata/database/gatk_bundle/b37/dbsnp_138.b37.vcf";
my $gatk_Mills_and_1000G_gold_standard_indels="/biodata/database/gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";
my $gatk_hapmap="/biodata/database/gatk_bundle/b37/hapmap_3.3.b37.vcf";
my $gatk_omni2="/biodata/database/gatk_bundle/b37/1000G_omni2.5.b37.vcf";
my $gatk_1000G_snp="/biodata/database/gatk_bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $gatk_b37_dbsnp_138="/biodata/database/gatk_bundle/b37/dbsnp_138.b37.vcf";
my $gatk_GetPileupSummaries=	"/biodata/database/gatk_bundle/b37/Mutect2/GetPileupSummaries/small_exac_common_3_b37.vcf.gz";
my $mutect_gnomad	=	"/biodata/database/gatk_bundle/b37/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz";



#global var
#my $cfg_dir = "$out_dir/cfg";	`mkdir -p $cfg_dir`;
my ( $bed_region, $gatk_interval_list, $bait_gatk_interval_list, $gatk_cnv_interval_list );
my $pon4gatk_vcf = "";

#常染色体Max_odds,常染色体Min_odds,X染色体Max_odds,X染色体Min_odds,Y染色体Max_odds,Y染色体Min_odds
my $bin_cnvkits		=	"$docker_cmd cnvkit.py";
my $bin_filterCNVkit	=	"$docker_cmd perl $bin_dir/cnvkit/filter_cnvkit.pl";
my $ref_cnvkits_access	=	"/biodata/software/CNVkit/cnvkit-master/data/access-5k-mappable.grch37.bed";
my $cnvkit_odds_limit_female =	"0.3,-0.4,0.3,-0.4,0.8,-3";
my $cnvkit_odds_limit_male =	"0.3,-0.4,0.8,-3,0.8,-3";
my $cnvkit_odds_limit_common =	$cnvkit_odds_limit_female;

my %h_germline_samples = ();	#存储所有call了germline变异的样本ID，最终综合统计一下性别等信息
my %h_bed_info = ("count"=>0);	#存储各bed数据结果
#所有normal结果存储
my %h_normal = ();
#读入数据配置文件，脚本打印


open IN, $infile or die $!;
my @lines = <IN>;
foreach my $line( @lines ){
	chomp $line;
	my @arr =split(/\t/,$line);
	my $t_modes = shift @arr;
	if( $t_modes =~ /target_bed/i ){#1. PE测序原始数据处理
		$gatk_interval_list = $arr[1];
		$bait_gatk_interval_list = $gatk_interval_list;
		&init_inFile($gatk_interval_list,$bait_gatk_interval_list);
	}
	elsif( $t_modes =~ /sample_pe/i ){#1. PE测序原始数据处理
		my ( $gender, $sample_id, $hpo_ids, $fq1, $fq2 ) = @arr;
		&init_inFile($fq1,$fq2);

		my ( $recal_bam,$sorted_bam,$table_getpileupsummaries,$markDup_metric ) = split(/,/,&gatk_preprocess( $sample_id,$fq1,$fq2));
		my ( $gatk_qc_prefix ) = split(/,/,&gatk_qc( $sample_id,$sorted_bam ));
		my ( $germline_raw_gvcf ) = split(/,/,&gatk4_germline_single_gvcf($sample_id,$recal_bam));
		my ( $germline_raw_vcf ) = split(/,/,&gatk4_germline_single_vcf($sample_id,$recal_bam));
		my ( $hardfilter_vcf ) = split(/,/,&gatk4_germline_single_hardfilter($sample_id,$germline_raw_vcf));
		&germline_vcf_anno( $sample_id, $hardfilter_vcf, "hardfilter" );
		my ( $vqsrfilter_vcf ) = split(/,/,&gatk4_germline_single_vqsrfilter($sample_id,$germline_raw_vcf));
		#cnnfilter 走不通暂时不用
		#my ( $cnnfilter_vcf ) = split(/,/,&gatk4_germline_single_cnn1D($sample_id,$germline_raw_vcf));
		#&germline_vcf_anno( $sample_id, $cnnfilter_vcf, "cnnfilter" );



	}
}
close IN;



sub gatk_preprocess{
		#步骤gatk preprocess
		my ( $sample_id, $in_fq1,$in_fq2) = @_;
		my $task_id = "$sample_id\_gatk_pre";

		#set names
		my $data_dir="$docker_out_dir/gatk_preprocess/$sample_id";
		my $temp_dir = "$data_dir/temp";
		my $trimFq_prefix="$data_dir/$sample_id.trim";
		my $bwa_prefix="$data_dir/$sample_id";
		my $raw_sam = "$bwa_prefix\.sam";
		my $markDup_bam = "$bwa_prefix\.markDup.bam";
		my $markDup_metric = "$bwa_prefix\.markDup.metrics.txt";
		my $sorted_bam="$bwa_prefix\.markDup.sorted.bam";
		my $recal_table="$bwa_prefix\.markDup.sorted.recal_data.table";
		my $recal_bam="$bwa_prefix\.markDup.sorted.recal.bam";
		my $table_getpileupsummaries="$recal_bam\.getpileupsummaries.table";


		&task_run(1,	$task_id, "$in_fq1,$in_fq2", "$recal_bam,$sorted_bam,$table_getpileupsummaries,$markDup_metric",
				"$bin_mkdir $data_dir",
				"$bin_mkdir $temp_dir",
				"$bin_bwa mem -t 16 -M $ref_fa $in_fq1 $in_fq2 >$raw_sam",
				"$bin_picard MarkDuplicates I=$raw_sam O=$markDup_bam M=$markDup_metric VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname",
				"$bin_picard AddOrReplaceReadGroups I=$markDup_bam O=$sorted_bam SORT_ORDER=coordinate RGID=$sample_id  RGLB=lb_$sample_id RGPL=illumina RGPU=pu_$sample_id RGSM=$sample_id TMP_DIR=$temp_dir",
				"$bin_samtools index $sorted_bam",
				"$bin_gatk BaseRecalibrator -R $ref_fa -I $sorted_bam -O $recal_table --known-sites $gatk_dbsnp --known-sites $gatk_Mills_and_1000G_gold_standard_indels --use-original-qualities",
				"$bin_gatk ApplyBQSR -R $ref_fa -I $sorted_bam -bqsr $recal_table -O $recal_bam  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30  --add-output-sam-program-record --create-output-bam-md5  --use-original-qualities",
				"$bin_gatk GetPileupSummaries -R $ref_fa -V $gatk_GetPileupSummaries --interval-set-rule INTERSECTION -L $gatk_interval_list -L $gatk_GetPileupSummaries -I $recal_bam -O $table_getpileupsummaries"
		);
}




#步骤gatk germline single sample call vcf
sub gatk4_germline_single_vcf{

		my ( $sample_id, $recal_bam ) = @_;
		my $task_id = "$sample_id\_CallVariantGermlineGATK4_single_vcf";
		my $data_dir="$docker_out_dir/germline_gatk4_single_vcf/$sample_id";
		#set output names
		my $germline_raw_vcf="$data_dir/germline.$sample_id\.raw.vcf.gz";
		my $germline_left_split_vcf="$data_dir/germline.$sample_id\.raw.left.split.vcf.gz";
		my $t_interval = "";
		$t_interval = "-L $gatk_interval_list" if( defined($gatk_interval_list) && $gatk_interval_list =~ /\S/ );

		#Hard filter following:https://software.broadinstitute.org/gatk/documentation/article?id=11069
		&task_run(1,	$task_id, "$recal_bam", "$germline_raw_vcf,$germline_left_split_vcf",
				"$bin_mkdir $data_dir",
				"$bin_gatk HaplotypeCaller -R $ref_fa -I $recal_bam $t_interval -O $germline_raw_vcf"
		);

}

#步骤gatk germline single sample call vcf
sub gatk4_germline_single_gvcf{

		my ( $sample_id, $recal_bam ) = @_;
		my $task_id = "$sample_id\_CallVariantGermlineGATK4_single_gvcf";
		my $data_dir="$docker_out_dir/germline_gatk4_single_gvcf/$sample_id";
		#set output names
		my $germline_raw_gvcf="$data_dir/germline.$sample_id\.raw.gvcf.gz";
		my $t_interval = "";
		$t_interval = "-L $gatk_interval_list" if( defined($gatk_interval_list) && $gatk_interval_list =~ /\S/ );

		#Hard filter following:https://software.broadinstitute.org/gatk/documentation/article?id=11069
		&task_run(1,	$task_id, "$recal_bam", "$germline_raw_gvcf",
				"$bin_mkdir $data_dir",
				"$bin_gatk HaplotypeCaller -R $ref_fa -I $recal_bam $t_interval -O $germline_raw_gvcf -ERC GVCF -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90"
		);

}


#步骤 gatk germline single sample, filter by hard filter
sub gatk4_germline_single_hardfilter{

		my ( $sample_id, $germline_raw_vcf  ) = @_;
		my $task_id = "$sample_id\_CallVariantGermlineGATK4_single_hardfilter";
		my $data_dir="$docker_out_dir/germline_gatk4_single_hardfilter/$sample_id";
		#set output names
		my $germline_filter_vcf="$data_dir/germline.$sample_id\.filter.vcf";
		my $germline_exclude_vcf="$data_dir/germline.$sample_id\.filter.exclude-non-variants.vcf";

		#Hard filter following:https://software.broadinstitute.org/gatk/documentation/article?id=11069
		&task_run(1,	$task_id, "$germline_raw_vcf", "$germline_exclude_vcf",
				"$bin_mkdir $data_dir",
				"$bin_gatk VariantFiltration -R $ref_fa -V $germline_raw_vcf -O $germline_filter_vcf --filter-name \\\"GATKHardFilter\\\" --filter-expression \\\"(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 3.0)\\\"",
				"$bin_gatk SelectVariants -R $ref_fa -V $germline_filter_vcf -O $germline_exclude_vcf --exclude-non-variants"
		);

}

#步骤gatk germline single sample
sub gatk4_germline_single_vqsrfilter{
		
		my ( $sample_id, $germline_raw_vcf  ) = @_;
		my $task_id = "$sample_id\_CallVariantGermlineGATK4_single_vqsrfilter";
		my $data_dir="$docker_out_dir/germline_gatk4_single_vqsrfilter/$sample_id";
		#my $germline_raw_gvcf="$data_dir/germline.$sample_id\.raw.gvcf.gz";
		my $germline_vqsr_pre="$data_dir/germline.$sample_id\.vqsr";
		my $germline_vr_snp_recal ="$germline_vqsr_pre.snp.recal";
		my $germline_vr_snp_tranches ="$germline_vqsr_pre.snp.tranches";
		my $germline_vr_snp_plot ="$germline_vqsr_pre.snp.plot.R";
		my $germline_vr_snp_VQSR ="$germline_vqsr_pre.snp.vcf";
		my $germline_vr_indel_recal ="$germline_vqsr_pre.indel.recal";
		my $germline_vr_indel_tranches ="$germline_vqsr_pre.indel.tranches";
		my $germline_vr_indel_plot ="$germline_vqsr_pre.indel.plot.R";
		my $germline_vr_indel_VQSR ="$germline_vqsr_pre.vcf";
		my $t_interval = "";
		$t_interval = "-L $gatk_interval_list" if( defined($gatk_interval_list) && $gatk_interval_list =~ /\S/ );
	
		&task_run(1,	$task_id, "$germline_raw_vcf", "$germline_vr_indel_VQSR",
				"$bin_mkdir $data_dir",
				"$bin_gatk VariantRecalibrator -R $ref_fa -V $germline_raw_vcf --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $gatk_hapmap --resource:omni,known=false,training=true,truth=true,prior=12.0 $gatk_omni2 --resource:1000G,known=false,training=true,truth=false,prior=10.0 $gatk_1000G_snp --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $gatk_b37_dbsnp_138     -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O $germline_vr_snp_recal --tranches-file $germline_vr_snp_tranches --rscript-file $germline_vr_snp_plot --max-gaussians 6",
				"$bin_gatk ApplyVQSR  -R $ref_fa -V $germline_raw_vcf -ts-filter-level 99.7 -mode SNP --recal-file $germline_vr_snp_recal --tranches-file $germline_vr_snp_tranches -O $germline_vr_snp_VQSR",
				"$bin_gatk VariantRecalibrator -R $ref_fa -V $germline_vr_snp_VQSR --resource:hapmap,known=false,training=true,truth=true,prior=12.0 $gatk_Mills_and_1000G_gold_standard_indels --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $gatk_b37_dbsnp_138 -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL -O $germline_vr_indel_recal --tranches-file $germline_vr_indel_tranches --rscript-file $germline_vr_indel_plot --max-gaussians 4",
				"$bin_gatk ApplyVQSR  -R $ref_fa -V $germline_vr_snp_VQSR -ts-filter-level 99.7 -mode INDEL --recal-file $germline_vr_indel_recal --tranches-file $germline_vr_indel_tranches -O $germline_vr_indel_VQSR"
		);
}


#步骤 gatk germline single sample, filter by cnn
#运行不通，报错
sub gatk4_germline_single_cnn1D{
		my ( $sample_id, $germline_raw_vcf  ) = @_;
		my $task_id = "$sample_id\_CallVariantGermlineGATK4_single_cnn1D";
		my $data_dir="$docker_out_dir/germline_gatk4_single_cnn1D/$sample_id";
		#set output names
		my $germline_anno_vcf="$data_dir/germline.$sample_id\.anno_cnn.vcf";
		my $germline_filter_vcf="$data_dir/germline.$sample_id\.cnn_filter.vcf";
		my $germline_exclude_vcf="$data_dir/germline.$sample_id\.filter.exclude-non-variants.vcf";
		my $germline_fix_vcf="$data_dir/germline.$sample_id\.filter.exclude-non-variants.vcf";
		my $snpeff_anno_stat="$data_dir/snpeff.$sample_id.stat.html";
		my $snpeff_anno_vcf="$data_dir/snpeff.$sample_id.vcf";
		my $snpeff_anno_fix_vcf="$data_dir/snpeff.$sample_id.fix.select_nm.annogene.vcf";
		my $annovar_anno_out="$data_dir/snpeff.$sample_id.fix.select_nm.annogene.hg19_anno";
		my $t_interval = "";
		$t_interval = "-L $gatk_interval_list" if( defined($gatk_interval_list) && $gatk_interval_list =~ /\S/ );

		#Hard filter following:
		&task_run(1,	$task_id, "$germline_raw_vcf", "$germline_filter_vcf",
				"$bin_mkdir $data_dir",
				"$bin_gatk CNNScoreVariants -R $ref_fa -V $germline_raw_vcf -O $germline_anno_vcf",
				"$bin_gatk FilterVariantTranches -V $germline_anno_vcf -O $germline_filter_vcf  --resource $gatk_hapmap  --resource $gatk_Mills_and_1000G_gold_standard_indels  --info-key CNN_1D    --snp-tranche 99.95  --indel-tranche 99.4");

}

#步骤 gatk germline single sample, filter by hard filter
sub germline_vcf_anno{

		my ( $sample_id, $germline_vcf, $tag  ) = @_;
		my $task_id = "$sample_id\_germline_vcf_anno_$tag";
		my $data_dir="$docker_out_dir/germline_vcf_anno_$tag/$sample_id";
		#set output names
		my $germline_lef_split_vcf="$data_dir/germline.$sample_id\.left.split.vcf";
		my $germline_fix_vcf="$data_dir/germline.$sample_id\.fix.vcf";
		my $snpeff_anno_stat="$data_dir/germline.$sample_id.fix.snpeff.stat.html";
		my $snpeff_anno_vcf="$data_dir/germline.$sample_id.fix.snpeff.vcf";
		my $snpeff_anno_fix_vcf="$data_dir/germline.$sample_id.fix.snpeff.select_nm.annogene.vcf";
		my $annovar_anno_out="$data_dir/germline.$sample_id.fix.snpeff.select_nm.annogene.hg19_annovar";

		#Hard filter following:https://software.broadinstitute.org/gatk/documentation/article?id=11069
		&task_run(1,	$task_id, "$germline_vcf", "$snpeff_anno_fix_vcf",
				"$bin_mkdir $data_dir",
				"$bin_bcftools norm -f $ref_fa -o $germline_lef_split_vcf -O v -m ->both -w 2000 $germline_vcf",
				"$bin_fixVcf4 -ref $ref_fa -in $germline_lef_split_vcf -out $germline_fix_vcf",
				"$bin_snpeff -s $snpeff_anno_stat $germline_fix_vcf >$snpeff_anno_vcf",
				"$docker_cmd cat $snpeff_anno_vcf|$bin_fix_snpeff|$bin_select_snpeff|$bin_anno_gene -clinvar_all $database_clinvar_stat_all -clinvar_mul $database_clinvar_stat_mul -hpo_g2p $database_g2p_hpo  -hpo_p2g $database_p2g_hpo >$snpeff_anno_fix_vcf",
				"$bin_table_annovar -out $annovar_anno_out -buildver hg19 -protocol clinvar_20190305,cosmic70,civic_20200317,civicRegion4SmallMut_20200317,gwasCatalog,avsnp150,gnomad211_genome,gnomad211_exome,dbnsfp35c -operation f,f,f,r,r,f,f,f,f $snpeff_anno_fix_vcf $annovar_humandb --remove --otherinfo --nastring . --vcfinput"
		);

}
#步骤gatk qc
sub gatk_qc{
		
		my ( $sample_id, $sorted_bam, $tag  ) = @_;
		my $task_id = "$sample_id\_gatk_qc";
		my $data_dir="$docker_out_dir/gatk_qc/$sample_id";
		#set output names
		my $qc_prefix="$data_dir/$sample_id.qc";
		#QC by picard & bedtools
		my $t_interval1 = "";
		my $t_interval2 = "";
		if( defined($gatk_interval_list) && $gatk_interval_list =~ /\S/ ){
			$t_interval1 = "INTERVALS=$gatk_interval_list";
			$bait_gatk_interval_list = $gatk_interval_list if(!defined($bait_gatk_interval_list));
			$t_interval2 = "BAIT_INTERVALS=$bait_gatk_interval_list TARGET_INTERVALS=$gatk_interval_list";
		}
		my $program = "PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectGcBiasMetrics PROGRAM=CollectSequencingArtifactMetrics PROGRAM=CollectQualityYieldMetrics";

		&task_run(1,	$task_id, "$sorted_bam", "$qc_prefix",
				"$bin_mkdir $data_dir",
				"$bin_picard CollectHsMetrics I=$sorted_bam R=$ref_fa O=$qc_prefix.hs_metrics.txt $t_interval2 PER_TARGET_COVERAGE=$qc_prefix.hs_metrics.per_target.txt COVERAGE_CAP=50000",
				"$bin_picard CollectMultipleMetrics $program I=$sorted_bam R=$ref_fa O=$qc_prefix.multiple_metrics $t_interval1"
		);

}






sub usageWithColor {
	print color "green";#change the text color
	print <<USAGE;
Description:
        Used to ...
Usage:  
        perl $0 [run_type] [options]
	run_type:
		create_pon4gatk	: 
        common Options:
	     -help : reveal help info
	     -in  <str>  : input file path, format to be defined
	     -out  <str>  : output dir path
        Example:
             perl $0 -in input.list -out out_dir/data_analysis >to_run.sh
Author & Contact:
	Mingming Liu
Last updated:
        2021-02-16
USAGE
	print color "reset";#change back the text color
}

