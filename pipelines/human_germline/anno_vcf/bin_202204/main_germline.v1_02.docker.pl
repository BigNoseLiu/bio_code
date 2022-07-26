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

my $print_flag = 1;
my ( $help, $infile, $out_dir, $GATK_PON, $optional_cmd, $debug );
#获取当前日期
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
$year+=1900;
$mon+=101;$mon =~ s/^1//;
$mday+=100;$mday =~ s/^1//;
$hour+=100;$hour=~ s/^1//;
$min+=100;$min=~ s/^1//;
$sec+=100;$sec=~ s/^1//;

#初始化输出文件路径名
#$out_dir = "data_$year$mon$mday$hour$min$sec";
my $inDir_toMount = "/in_dir_mount";
my %h_cmd_log = ();		#存储各样本命令
my %h_sample_info = ();		#存储各样本数据结果
my %h_infile4docker = ();
my %h_infile4preTask = ();	#存储文件的生成任务


GetOptions(
	"help"=>\$help,
	"debug"=>\$debug,
	"opt"=>\$optional_cmd,
	"in=s"=>\$infile,
	"out=s"=>\$out_dir,
	"pon=s"=>\$GATK_PON,
);

#打印帮助文件
if ( defined($help) || !defined($infile) || !defined($out_dir) ) {
	&usageWithColor();
	exit 0;
}

#加前缀路径，并创建文件夹

my $run_type = $ARGV[0];
sub err_print{
	print STDERR color "red";#change the text color
	foreach my $t_cmd( @_ ){
		$t_cmd = &linkFileToDocker($t_cmd);
		print STDERR $t_cmd."\n";
	}
	print STDERR color "reset";#change back the text color
}

#print out cmd lines
sub cmd_print{
	#my $print_flag = shift@_;
	my $print_flag = 1;
	foreach my $t_cmd( @_ ){
		$t_cmd = &linkFileToDocker($t_cmd);
		if( !defined( $debug ) ){
			$t_cmd =~ s/^\s*(\S)/$1/;
			$t_cmd =~ s/(\S)\s*$/$1/;
			$t_cmd =~ s/(biojupyter)/$1 bash -c "/;
			$t_cmd =~ s/$/"/;
			if( $print_flag == 1 ){
				print $t_cmd."\n";
				if( $t_cmd =~ /^#&J/ ){
				#	print 'echo -n "'.$t_cmd.'	;'.'date "+%Y-%m-%d %H:%M:%S";'."\t";#打印到STDOUT
				#	print 'echo -n "'.$t_cmd.'	>&2;'.'date "+%Y-%m-%d %H:%M:%S" >&2'."\n";#打印到STDERR
				}
			}
		}
	}
}


#处理文件到docker可识别的方式
sub linkFileToDocker{
	my $cmd = shift @_;
	$cmd = " $cmd ";
	foreach my $raw_file( sort {$a cmp $b} keys(%h_infile4docker) ){
			#输入文件路径改为dockers识别的路径
			if(  $cmd =~ /docker\s+run\s/ && $cmd =~ /([=\s])$raw_file([=\s])/ ){
				my $in_file = $raw_file;	$in_file =~ s/^\//$inDir_toMount\//;
				$cmd =~ s/([=\s])$raw_file([\s=])/$1$in_file$2/g;
				$cmd =~ s/docker(\s+)run/docker run -v $raw_file:$in_file /;
			}
	}
	return $cmd;
}

sub init_inFile{
	foreach my $file_path(@_){
		$h_infile4preTask{$file_path} = "NA";
		$h_infile4docker{$file_path}++;
	}
}

#moreover, optional print out lines
sub opt_print{
	foreach my $t_cmd( @_ ){
		$t_cmd = &linkFileToDocker($t_cmd);
		if( defined( $optional_cmd ) ){
			print $t_cmd."\n";
		}
	}
}

sub task_admin{
	
	#初始化
	my ($task_id, $sample_id, $in_files, $out_files)= @_;
	my %pre_ids = ("NA"=>1);
	my @arr_return = (0);


	#判断该任务是否已运行过
	if( defined( $h_cmd_log{$task_id} ) ){
		$arr_return[0] = 1;
		return @arr_return;
	}
	$h_cmd_log{$task_id}++;

	#输入文件处理
	if( $in_files ne "-" ){
		foreach my $file_type( split(/,/,$in_files) ){
			my $file_path = $h_sample_info{"sample"}{$sample_id}{$file_type}{"path"};
			push @arr_return,$file_path;
			$pre_ids{ $h_infile4preTask{$file_path} } ++;
		}
	}
	my $str_pre_id = join(",",sort {$a cmp $b} keys(%pre_ids));
	$str_pre_id =~ s/,NA,//;	$str_pre_id =~ s/^NA,//;	$str_pre_id =~ s/,NA$//;
	push @arr_return, $str_pre_id;

	#输出文件处理
	if( $out_files ne "-" ){
		foreach my $name2file( split(/,/,$out_files) ){
			if( $name2file =~ /^\s*(\S+):(\S+)\s*$/ ){
				my ($file_type,$file_path) = ($1,$2);
				$h_sample_info{"sample"}{$sample_id}{$file_type}{"path"} = $file_path;
				$h_infile4preTask{$file_path} = $task_id;
			}
			else{
				print STDERR "Illegal out_file format $name2file\n";
				exit(0);
			}
		}
	}


	return @arr_return;#样例：my ($flag_dup, $raw_fq1, $raw_fq2, $pre_task)  = &task_admin($task_id, $sample_id, "raw_fq1,raw_fq2", "clean_fq1:$out_trimFq1,clean_fq2:$out_trimFq2");
}




#only print out debug lines
sub debug_print{
	foreach my $t_cmd( @_ ){
		$t_cmd = &linkFileToDocker($t_cmd);
		if( defined( $debug ) ){
			print $t_cmd."\n";
		}
	}
}

my $bin_dir = '';#"/biodata/biopipe/human_germline/anno_vcf/bin_202204";
my $biodata_dir='';#'/data/backup/home/liumingming/biodata/';
if( $Bin =~ /^(.*)\/biodata\/(.*)$/ ){
	$bin_dir = "/biodata/$2";
	$biodata_dir="$1/biodata";
}
else{
	print STDERR "Err: you're not allowed to change the name of biodata!!!\n";
	exit(0);
}

#设置输出目录
`mkdir -p $out_dir $out_dir/temp`;
my $docker_out_dir = "/mnt/$out_dir";
my $temp_dir = "$docker_out_dir/temp";

#设置docker运行命令
my $ug_id='1020:1022';
my $docker_name='liumingming1988/biodocker';
my $docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u $ug_id $docker_name ";

sub docker_print{
	foreach my $t_cmd( @_ ){
		print $docker_cmd.' bash -c "'.$t_cmd.'"'."\n";
	}
}




my $bin_mkdir		=	"$docker_cmd mkdir -p";
#软件变量配置
my $bin_fastp		=	"$docker_cmd fastp";
my $bin_trimmomatic	=	"$docker_cmd trimmomatic";
my $bin_polysolver="/usr/bin/docker run -P --name \$NAME -v \$DIR1:\$DIR2 sachet/polysolver:v4 bash /home/polysolver/scripts/shell_call_hla_type \$DIR2/\$BAM Unknown 1 hg19 STDFQ 0 \$DIR2/\$OUT_DIR; /usr/bin/docker rm \$NAME";
my $bin_statHLA="perl $bin_dir/HLA/stat_HLA.v1.01.pl";
#my $bin_polysolver="/usr/bin/docker run -P --name \$NAME -v \$DIR1:\$DIR2 sachet/polysolver:v4 bash /home/polysolver/scripts/shell_call_hla_type \$DIR2/\$BAM Unknown 1 hg19 STDFQ 0 \$DIR2; /usr/bin/docker rm \$NAME";
my $bin_msisensor	=	"$docker_cmd msisensor-pro";
my $bin_factera		=	"$docker_cmd perl /biodata/software/FACTERA/v1.4.4/factera.pl";
my $factera_bed		=	"/biodata/software/FACTERA/v1.4.4/exons.noChr_noMT.bed";
my $factera_bed2	=	"/biodata/software/FACTERA/v1.4.4/exons.noChr_noMT.exonlevel.bed";
my $bin_facterastat	=	"perl $bin_dir/factera/stat_factera.pl";
my $bin_filter_umi_cfbestV2=	"perl $bin_dir/umi/trimfq_cfbest_umi_v2.pl";
my $bin_fgbio		=	"$docker_cmd fgbio";
my $bin_varscan		=	"$docker_cmd varscan";
my $bin_snpeff		=	"java -jar /biodata/git_code/software/snpeff/snpEff5_0/snpEff.jar -v GRCh37.p13.RefSeq";
my $bin_select_snpeff	=	"perl $bin_dir/snpeff/select_standard_hgvs_from_snpeff_Eff.v1.01.pl";
my $bin_fix_snpeff	=	"perl $bin_dir/snpeff/fix_snpeff_Eff.v1.pl";
my $bin_bedtools	=	"$docker_cmd bedtools";
my $bin_qualimap="$bin_dir/../../software/QualiMap/qualimap_v2.2.1/qualimap";
my $bin_filterMutect	=	"$docker_cmd perl $bin_dir/gatk4_somatic_SNVIndel/filter_mutect2_pair.leftSplit.pl";
my $bin_fixVcf4		=	"$docker_cmd perl $bin_dir/gatk4_germline_SNVIndel/fix_gatk4_vcf4_2.pl";
my $bin_tag_nearbyOverlap_variant="perl $bin_dir/gatk4_germline_SNVIndel/tag_nearbyOverlap_variant.v1.01.pl";
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
my $bin_exomiser	=	"java -Xms2g -Xmx4g -jar /biodata/git_code/software/exomiser/exomiser-cli-13.0.0/exomiser-cli-13.0.0.jar";
my $config_exomiser	=	"/biodata/git_code/software/exomiser/exomiser-cli-13.0.0/application.properties";
my $bin_create_exomiser_exome=	"cat $bin_dir/exomiser/exome.yml|perl $bin_dir/exomiser/create_yml.pl";
#annovar
my $bin_annovar		=	"/biodata/git_code/software/annovar/bin";
my $annovar_humandb	=	"/biodata/databases/annovar_humandb";
#my $annovar_para	= 	'-protocol clinvar_20210501,HGMD2017,gwasCatalog,avsnp150,gnomad211_genome,gnomad211_exome,esp6500siv2_all,dbnsfp42a,dbscsnv11,1000g2014oct_all,1000g2014oct_eas,1000g2014oct_eur,1000g2014oct_amr,1000g2014oct_afr,1000g2014oct_sas -operation f,f,r,f,f,f,f,f,f,f,f,f,f,f,f';
my $annovar_para	= 	'-protocol clinvar_20220416,clinvar_samepos_20220416,clinvar_flank5_20220416,HGMD2017_samepos,HGMD2017_flank5,avsnp150,gnomad211_genome,gnomad211_exome,exac03,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_afr,1000g2015aug_sas,dbnsfp42a,dbscsnv11 -operation f,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f';
my $bin_annovar_withPara=	"perl $bin_annovar/table_annovar.pl -buildver hg19 $annovar_para --remove --otherinfo --nastring . --vcfinput";
#my $bin_annovar_withPara=	"perl $bin_annovar/table_annovar.pl -buildver hg19 -protocol refGeneWithVer,cosmic70,cosmic81,civic_20200317,civicRegion4SmallMut_20200317,clinvar_20190305,HGMD2017,gwasCatalog,avsnp150,gnomad211_genome,gnomad211_exome,esp6500siv2_all,dbnsfp35a,nci60,1000g2014oct_all,1000g2014oct_eas,1000g2014oct_eur,1000g2014oct_amr,1000g2014oct_afr,1000g2014oct_sas,ljb26_all -operation gx,f,f,f,r,f,f,r,f,f,f,f,f,f,f,f,f,f,f,f,f --remove --otherinfo --nastring . --vcfinput";
my $bin_anno_gene	=	"perl $bin_dir/var_inter/anno_gene.v2.01.pl";
#HPO:gene to disease
my $database_g2p_hpo	=	"/biodata/databases/HPO/20200326/genes_to_phenotype.txt";
my $database_p2g_hpo	=	"/biodata/databases/HPO/20200326/phenotype_to_genes.txt";
#ClinVar, note: must be matched with $bin_snpeff(eg. -v GRCh37.p13.RefSeq)
my $database_clinvar_stat_all=	"/biodata/databases/clinvar/20200326/clinvar_20200316.snpeff.fix.vcf.filter_no.stat.xls";
my $database_clinvar_stat_mul=	"/biodata/databases/clinvar/20200326/clinvar_20200316.snpeff.fix.vcf.filter_multiple.stat.xls";
#ref genome
my $ref_fa		=	"/biodata/databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta";
#my $ref_fa		=	"/biodata/database/genomes/b37_decoy/human_g1k_v37_decoy.fasta";
my $ref_fa_dict		=	"/biodata/databases/genomes/b37_decoy/human_g1k_v37_decoy.dict";
my $ref_fa_2bit		=	"/biodata/databases/genomes/b37_decoy/human_g1k_v37_decoy.2bit";
my $ref_flat		=	"/biodata/software/CNVkit/database/hg19/GRch37.refFlat.txt";
my $ref_msi_microsatellite=	"/biodata/databases/genomes/b37_decoy/human_g1k_v37_decoy.microsatellites.list";
#varDict
my $bin_vardict="$bin_dir/../../software/VarDict/VarDictJava-master/VarDict/VarDict-1.8.0/bin/VarDict";
my $bin_vardict_teststrandbias="$bin_dir/../../software/VarDict/VarDict-master/teststrandbias.R";
my $bin_vardict_var2vcf_valid="perl $bin_dir/../../software/VarDict/VarDict-master/var2vcf_valid.pl";
#gatk
my $bin_gatk		=	"gatk";
my $gatk_dbsnp		=	"/biodata/databases/gatk_bundle/b37/dbsnp_138.b37.vcf";
my $gatk_hapmap="/biodata/databases/gatk_bundle/b37/hapmap_3.3.b37.vcf";
my $gatk_omni2="/biodata/databases/gatk_bundle/b37/1000G_omni2.5.b37.vcf";
my $gatk_1000G_snp="/biodata/databases/gatk_bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $gatk_b37_dbsnp_138="/biodata/databases/gatk_bundle/b37/dbsnp_138.b37.vcf";
my $gatk_Mills_and_1000G_gold_standard_indels=	"/biodata/databases/gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";
my $gatk_b37_mills			=	"/biodata/databases/gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";
my $gatk_GetPileupSummaries=	"/biodata/databases/gatk_bundle/b37/Mutect2/GetPileupSummaries/small_exac_common_3_b37.vcf.gz";
my $mutect_gnomad	=	"/biodata/databases/gatk_bundle/b37/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz";



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
chomp @lines;
#&germline_pip(@lines);
foreach my $line( @lines ){
	next if($line =~ /^#/);
	my @arr =split(/\t/,$line);
	my ( $sample_id,$hpo_ids, $to_anno_vcf ) = @arr;
	$to_anno_vcf = "/mnt/$to_anno_vcf";
	my $anno_type = "germline_anno";
	my $data_dir="$docker_out_dir/$anno_type/$sample_id";
	#print STDERR "$anno_type\t$sample_id\t$to_anno_vcf\n";
	my $left_align_vcf = "$data_dir/$sample_id.final.leftalign.vcf";
	my $yml_exomiser_file = "$data_dir/exomiser.$sample_id.yml";
	my $exomiser_anno_prefix ="$data_dir/exomiser.$sample_id";
	my $snpeff_anno_vcf="$data_dir/snpeff.$sample_id.vcf";
	my $snpeff_anno_stat="$data_dir/snpeff.$sample_id.stat.html";
	my $anno_gene_vcf="$data_dir/$sample_id.final.leftalign.annoGene.vcf";
	my $annoVar_out="$data_dir/$sample_id.final.leftalign.annoGene.annoVar";

	my $last_job = "$sample_id\_CallVariant$anno_type";
	&cmd_print("#&J	$sample_id\_AnnoVariant_$anno_type	$last_job");
	my $anno_job_ids .= ",$sample_id\_AnnoVariant_$anno_type";
	&docker_print("mkdir -p $data_dir");
	&docker_print("$bin_gatk LeftAlignAndTrimVariants --split-multi-allelics -R $ref_fa -V $to_anno_vcf -O $left_align_vcf");
	my $anno_exomiser_cmd = "";
	if( $hpo_ids ne "-" ){
		&docker_print("$bin_create_exomiser_exome liumm_vcf_file_path:$left_align_vcf liumm_output_prefix:$exomiser_anno_prefix liumm_hpo_ids:$hpo_ids >$yml_exomiser_file");
		&docker_print("$bin_exomiser --analysis $yml_exomiser_file --spring.config.location=$config_exomiser");
		$anno_exomiser_cmd = "-exomiser $exomiser_anno_prefix.vcf";
	}
	&docker_print("$bin_snpeff -s $snpeff_anno_stat $left_align_vcf >$snpeff_anno_vcf");
	&docker_print("cat $snpeff_anno_vcf|$bin_fix_snpeff|$bin_select_snpeff|$bin_tag_nearbyOverlap_variant -ref $ref_fa >$anno_gene_vcf");#|$bin_anno_gene $anno_exomiser_cmd -clinvar_all $database_clinvar_stat_all -clinvar_mul $database_clinvar_stat_mul -hpo_g2p $database_g2p_hpo -hpo_p2g $database_p2g_hpo >$anno_gene_vcf");
	#注释各类数据库信息by annovar
	&docker_print("$bin_annovar_withPara -out $annoVar_out $anno_gene_vcf $annovar_humandb");
	#$vcf_files_to_stat .= " $sample_id\_$anno_type:$anno_type:$annoVar_out\.hg19_multianno.vcf";
}
close IN;

&QC_stat();
&gatk_join_call();


#分析Germline变异流程
sub germline_pip{
	foreach my $line( @_ ){
		chomp $line;
		my @arr =split(/\t/,$line);
		my $t_modes = shift@arr;

		#获取目标区域bed文件
		if( $t_modes =~ /^#target_bed\s*$/i ){
			$bed_region = $arr[0];
			$gatk_interval_list = $arr[1];
			&init_inFile($bed_region, $gatk_interval_list);
		}
		elsif( $t_modes =~ /^#bait_bed\s*$/i ){
			$bait_gatk_interval_list = $arr[1];
			&init_inFile($bait_gatk_interval_list);
		}
		elsif( $t_modes =~ /sample_pe/i ){#1. PE测序原始数据处理
			my ( $gender, $sample, $hpo_ids, $fq1, $fq2 ) = @arr;
			if( defined($h_sample_info{"sample"}{$sample}) ){
				&err_print("Err: duplicate sample $sample");
				exit(0);
			}

			$h_sample_info{"sample"}{$sample}{"id"} = $sample;
			$h_sample_info{"sample"}{$sample}{"gender"} = $gender;
			$h_sample_info{"sample"}{$sample}{"raw_fq1"}{"path"} = $fq1;
			$h_sample_info{"sample"}{$sample}{"raw_fq2"}{"path"} = $fq2;
			$h_sample_info{"sample"}{$sample}{"clean_fq1"}{"path"} = $fq1;
			$h_sample_info{"sample"}{$sample}{"clean_fq2"}{"path"} = $fq2;
			&init_inFile($fq1, $fq2);
			#&fastp_qc(0,$sample);
			&gatk_preprocess(1,$sample);
			&gatk_qc(1,$sample);
			&germline_gatk4(1,$sample);
		}
	}
}



#1. preprocess
sub fastp_qc{
	
	#任务ID
	my ($print_flag, $sample_id) = @_;
	my $task_id = "$sample_id\_fastp_QC";
	#初始化文件名
	my $data_dir="$out_dir/fastp/$sample_id";
	my $trimFq_prefix="$data_dir/$sample_id.trim";
	my $out_trimFq1 = "$trimFq_prefix.R1.fq.gz";
	my $out_trimFq2 = "$trimFq_prefix.R2.fq.gz";

	my ($flag_dup, $raw_fq1, $raw_fq2, $pre_task)  = &task_admin($task_id, $sample_id, "raw_fq1,raw_fq2", "clean_fq1:$out_trimFq1,clean_fq2:$out_trimFq2");
	return if($flag_dup == 1);
	
	#运行脚本
	&cmd_print($print_flag, "#&J	$task_id	$pre_task");
	&cmd_print($print_flag, "$bin_mkdir $data_dir ");
	&cmd_print($print_flag, "$bin_fastp -i $raw_fq1 -I $raw_fq2 -o $trimFq_prefix.R1.fq.gz -O $trimFq_prefix.R2.fq.gz  -h $trimFq_prefix.html -j $trimFq_prefix.json");

}


sub gatk_preprocess{

	#任务ID
	my ($print_flag, $sample_id) = @_;
	my $task_id = "$sample_id\_gatk_pre";

	#set names
	my $data_dir="$out_dir/gatk_preprocess/$sample_id";
	my $trimFq_prefix="$data_dir/$sample_id.trim";
	my $bwa_prefix="$data_dir/$sample_id";
	my $raw_sam = "$bwa_prefix\.sam";
	my $markDup_bam = "$bwa_prefix\.markDup.bam";
	my $markDup_metric = "$bwa_prefix\.markDup.metrics.txt";
	my $sorted_bam="$bwa_prefix\.markDup.sorted.bam";
	my $recal_table="$bwa_prefix\.markDup.sorted.recal_data.table";
	my $recal_bam="$bwa_prefix\.markDup.sorted.recal.bam";
	my $table_getpileupsummaries="$recal_bam\.getpileupsummaries.table";
	
	my ($flag_dup, $in_fq1, $in_fq2, $pre_task)  = &task_admin($task_id, $sample_id, "clean_fq1,clean_fq2", "recal_bam:$recal_bam,sorted_bam:$sorted_bam,table_getpileupsummaries:$table_getpileupsummaries,dup_metric:$markDup_metric");
	return if($flag_dup == 1);

	#preprocess following GATK best practice
	&cmd_print($print_flag, "#&J	$task_id	$pre_task\n");
	&cmd_print($print_flag, "$bin_mkdir $data_dir");
	&cmd_print($print_flag, "$bin_bwa mem -t 16 -M $ref_fa $in_fq1 $in_fq2 >$raw_sam");
	&cmd_print($print_flag, "$bin_picard MarkDuplicates I=$raw_sam O=$markDup_bam M=$markDup_metric VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname");
	&cmd_print($print_flag, "$bin_picard AddOrReplaceReadGroups I=$markDup_bam O=$sorted_bam SORT_ORDER=coordinate RGID=$sample_id  RGLB=lb_$sample_id RGPL=illumina RGPU=pu_$sample_id RGSM=$sample_id TMP_DIR=$temp_dir");
	&cmd_print($print_flag, "$bin_samtools index $sorted_bam");
	&cmd_print($print_flag, "$bin_gatk BaseRecalibrator -R $ref_fa -I $sorted_bam -O $recal_table --known-sites $gatk_dbsnp --known-sites $gatk_Mills_and_1000G_gold_standard_indels --use-original-qualities");
	&cmd_print($print_flag, "$bin_gatk ApplyBQSR -R $ref_fa -I $sorted_bam -bqsr $recal_table -O $recal_bam  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30  --add-output-sam-program-record --create-output-bam-md5  --use-original-qualities");
	&cmd_print($print_flag, "$bin_gatk GetPileupSummaries -R $ref_fa -V $gatk_GetPileupSummaries --interval-set-rule INTERSECTION -L $gatk_interval_list -L $gatk_GetPileupSummaries -I $recal_bam -O $table_getpileupsummaries");

}

sub gatk_qc{

	my ($print_flag, $sample_id) = @_;;
	my $task_id = "$sample_id\_gatk_qc";

	#set output names
	my $data_dir="$out_dir/gatk_qc/$sample_id";
	my $qc_prefix="$data_dir/$sample_id.picardQC";

	my ($flag_dup, $sorted_bam, $pre_task )  = &task_admin($task_id, $sample_id, "sorted_bam", "gatk_qc:$qc_prefix*metrics*");
	return if($flag_dup == 1);

	#QC by picard & bedtools
	my $t_interval1 = "";
	my $t_interval2 = "";
	if( defined($gatk_interval_list) && $gatk_interval_list =~ /\S/ ){
		$t_interval1 = "INTERVALS=$gatk_interval_list";
		$bait_gatk_interval_list = $gatk_interval_list if(!defined($bait_gatk_interval_list));
		$t_interval2 = "BAIT_INTERVALS=$bait_gatk_interval_list TARGET_INTERVALS=$gatk_interval_list";
	}
	$h_sample_info{"sample"}{$sample_id}{"QC"}{"$qc_prefix*metrics*"}++;
	my $program = "PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectGcBiasMetrics PROGRAM=CollectSequencingArtifactMetrics PROGRAM=CollectQualityYieldMetrics";

	#打印命令
	&cmd_print($print_flag, "#&J	$task_id	$pre_task\n");
	&cmd_print($print_flag, "$bin_mkdir $data_dir");
	#&cmd_print($print_flag, "$bin_picard CollectHsMetrics I=$sorted_bam R=$ref_fa O=$qc_prefix.hs_metrics.txt $t_interval2 PER_TARGET_COVERAGE=$qc_prefix.hs_metrics.per_target.txt COVERAGE_CAP=50000");
	&cmd_print($print_flag, "$bin_picard CollectHsMetrics MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 CLIP_OVERLAPPING_READS=false I=$sorted_bam R=$ref_fa O=$qc_prefix.hs_metrics.txt $t_interval2 PER_TARGET_COVERAGE=$qc_prefix.hs_metrics.per_target.txt COVERAGE_CAP=50000");
	#&cmd_print($print_flag, "$bin_picard CollectMultipleMetrics $program I=$sorted_bam R=$ref_fa O=$qc_prefix.multiple_metrics $t_interval1");
	#&cmd_print("sed '1,11d' $qc_prefix.hs_metrics.txt|awk '\$2>0' > $qc_prefix.hs_metrics.plot_data.txt");
	#&opt_print("$bin_picard CollectHsMetrics I=$sorted_bam R=$ref_fa O=$qc_prefix.hs_metrics.mq0_bq0.txt $t_interval2 PER_TARGET_COVERAGE=$qc_prefix.hs_metrics.per_target.mq0_bq0.txt COVERAGE_CAP=50000 MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0");
	#&opt_print("$bin_picard CollectWgsMetrics I=$sorted_bam R=$ref_fa O=$qc_prefix.wgs_metrics.txt $t_interval1 COVERAGE_CAP=50000 MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 COUNT_UNPAIRED=true");

	#&cmd_print("$bin_qualimap bamqc -bam $sorted_bam -gff $bed_region -outdir $data_dir -outformat pdf");
	#$h_sample_info{"sample"}{$sample_id}{"QC"} = $data_dir."/genome_results.txt";
}

sub QC_stat{

	my $task_id = "AllStatQC";

	my $pre_task = "";
	my $input_str = "";
	foreach my $sample_id( keys(%{$h_sample_info{"sample"}}) ){
		foreach my $metric_id('dup_metric','gatk_qc'){
			if( defined($h_sample_info{"sample"}{$sample_id}{$metric_id}) ){
				$pre_task .= ",".$h_infile4preTask{$h_sample_info{"sample"}{$sample_id}{$metric_id}{"path"}};
				$input_str .= " $sample_id:".$h_sample_info{"sample"}{$sample_id}{$metric_id}{"path"};
			}
		}
	}

	if($pre_task =~ s/^,//){
		my $data_dir = "$out_dir/final_result";
		&cmd_print($print_flag, "#&J	$task_id	$pre_task\n");
		&cmd_print($print_flag, "$bin_mkdir $data_dir");
		#&cmd_print("$bin_statQC $data_dir/QC_stat.xls $input_str");
		&cmd_print($print_flag, "$bin_statPicardQC $data_dir/AllStat.QC_stat.xls $input_str");
	}
}


#Germline变异检测by GATK4
sub germline_gatk4{

	my ($print_flag, $sample_id) = @_;
	my $task_id = "$sample_id\_CallVariantGermlineGATK4";

	#set output names
	my $data_dir="$out_dir/germline_gatk4/$sample_id";
	my $germline_raw_vcf="$data_dir/germline.$sample_id\.raw.vcf.gz";
	my $germline_raw_gvcf="$data_dir/germline.$sample_id\.raw.g.vcf.gz";
	my $germline_vqsr_pre="$data_dir/germline.$sample_id\.raw.vqsr";
	my $germline_vr_snp_recal ="$germline_vqsr_pre.snp.recal";
	my $germline_vr_snp_tranches ="$germline_vqsr_pre.snp.tranches";
	my $germline_vr_snp_plot ="$germline_vqsr_pre.snp.plot.R";
	my $germline_vr_snp_VQSR ="$germline_vqsr_pre.snp.vcf";
	my $germline_vr_indel_recal ="$germline_vqsr_pre.indel.recal";
	my $germline_vr_indel_tranches ="$germline_vqsr_pre.indel.tranches";
	my $germline_vr_indel_plot ="$germline_vqsr_pre.indel.plot.R";
	my $germline_vr_indel_VQSR ="$germline_vqsr_pre.vcf";
	my $germline_filter_vcf="$data_dir/germline.$sample_id\.filter.vcf";
	my $germline_exclude_vcf="$data_dir/germline.$sample_id\.filter.exclude-non-variants.vcf";
	
	$h_germline_samples{$sample_id}++;#存储所有call了germline变异的样本ID，最终综合统计一下性别等信息

	my ($flag_dup, $recal_bam, $pre_task )  = &task_admin($task_id, $sample_id, "recal_bam", "final_vcf:$germline_exclude_vcf,gvcf:$germline_raw_gvcf");
	return if($flag_dup == 1);
	
	#germline variant call with GATK4
	my $t_interval = "";
	$t_interval = "-L $gatk_interval_list" if( defined($gatk_interval_list) && $gatk_interval_list =~ /\S/ );

	&cmd_print($print_flag, "#&J	$task_id	$pre_task\n");
	&cmd_print($print_flag, "$bin_mkdir $data_dir");
	&cmd_print($print_flag, "$bin_gatk HaplotypeCaller -R $ref_fa -I $recal_bam $t_interval -O $germline_raw_gvcf -ERC GVCF -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90");
	&cmd_print($print_flag, "$bin_gatk HaplotypeCaller -R $ref_fa -I $recal_bam $t_interval -O $germline_raw_vcf");
	&cmd_print($print_flag, "$bin_gatk VariantRecalibrator -R $ref_fa -V $germline_raw_vcf --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $gatk_hapmap --resource:omni,known=false,training=true,truth=true,prior=12.0 $gatk_omni2 --resource:1000G,known=false,training=true,truth=false,prior=10.0 $gatk_1000G_snp --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $gatk_b37_dbsnp_138     -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O $germline_vr_snp_recal --tranches-file $germline_vr_snp_tranches --rscript-file $germline_vr_snp_plot --max-gaussians 6");
	&cmd_print($print_flag, "$bin_gatk ApplyVQSR  -R $ref_fa -V $germline_raw_vcf -ts-filter-level 99.7 -mode SNP --recal-file $germline_vr_snp_recal --tranches-file $germline_vr_snp_tranches -O $germline_vr_snp_VQSR");

	&cmd_print($print_flag, "$bin_gatk VariantRecalibrator -R $ref_fa -V $germline_vr_snp_VQSR --resource:hapmap,known=false,training=true,truth=true,prior=12.0 $gatk_b37_mills --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $gatk_b37_dbsnp_138 -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0              -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL -O $germline_vr_indel_recal --tranches-file $germline_vr_indel_tranches --rscript-file $germline_vr_indel_plot --max-gaussians 4");
	&cmd_print($print_flag, "$bin_gatk ApplyVQSR  -R $ref_fa -V $germline_vr_snp_VQSR -ts-filter-level 99.7 -mode INDEL --recal-file $germline_vr_indel_recal --tranches-file $germline_vr_indel_tranches -O $germline_vr_indel_VQSR");

	#&cmd_print($print_flag, "$bin_gatk VariantFiltration -R $ref_fa -V $germline_raw_vcf -O $germline_filter_vcf --filter-name \\\"GATKHardFilter\\\" --filter-expression \\\"(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 3.0)\\\"");#Hard filter following:https://software.broadinstitute.org/gatk/documentation/article?id=11069
	#&cmd_print($print_flag, "$bin_gatk SelectVariants -R $ref_fa -V $germline_filter_vcf -O $germline_exclude_vcf --exclude-non-variants");

	#&cmd_print("$bin_fixVcf4 -ref $ref_fa -in $germline_filter_vcf.exclude-non-variants.vcf -out $germline_filter_vcf.exclude-non-variants.fix.vcf");
	#left align (not necessary for single sample)
	#&cmd_print("$bin_bcftools norm -f $ref_fa $germline_filter_vcf.exclude-non-variants.fix.vcf >$germline_filter_vcf.exclude-non-variants.fix.leftalign.vcf --multiallelics -");
	#&cmd_print("$bin_bcftools norm -f $ref_fa $germline_filter_vcf.exclude-non-variants.fix.leftalign.vcf >$final_vcf -d none");

}
sub gatk_join_call{
	my %h_gvcfs = ();
	my $gvcfs = "";
	foreach my $sample_id( keys(%{$h_sample_info{"sample"}})){
		if( defined($h_sample_info{"sample"}{$sample_id}{"gvcf"}) ){
			$h_gvcfs{$h_sample_info{"sample"}{$sample_id}{"gvcf"}{"path"}}++;
			$gvcfs .= " --INPUT ".$h_sample_info{"sample"}{$sample_id}{"gvcf"}{"path"};
		}
	}
	if( $gvcfs =~ /\S/ ){
		my $data_dir="$out_dir/germline_gatk4/join.g.vcf.gz";
		&cmd_print($print_flag, "$bin_gatk MergeVcfs $gvcfs --OUTPUT $data_dir");
	}
}

#Small Variant变异注释
sub anno_germline_variant{

	my ($print_flag, $sample_id) = @_;
	my $task_id = "$sample_id\_anno_germline";

	#set names
	my $data_dir="$out_dir/germline_anno";
	my $yml_exomiser_file = "$data_dir/exomiser.$sample_id.yml";	$h_sample_info{"sample"}{$sample_id}{"yml_exomiser_file"} = $yml_exomiser_file;
	my $exomiser_anno_prefix ="$data_dir/exomiser.$sample_id";	$h_sample_info{"sample"}{$sample_id}{"exomiser_anno_tsv"} = $exomiser_anno_prefix.".variants.tsv";
	my $snpeff_anno_vcf="$data_dir/snpeff.$sample_id.vcf";	$h_sample_info{"sample"}{$sample_id}{"snpeff_anno_vcf"} = $snpeff_anno_vcf;
	my $snpeff_anno_stat="$data_dir/snpeff.$sample_id.stat.html";	$h_sample_info{"sample"}{$sample_id}{"snpeff_anno_stat"} = $snpeff_anno_stat;
	my $annovar_anno_out ="$data_dir/annovar.$sample_id";	$h_sample_info{"sample"}{$sample_id}{"annovar_anno_out"} = $annovar_anno_out;
	my $final_anno = "$data_dir/final.$sample_id.annoResult.txt";	$h_sample_info{"sample"}{$sample_id}{"final_anno"} = $final_anno;


	my ($flag_dup, $final_vcf, $pre_task)  = &task_admin($task_id, $sample_id, "final_vcf", "yml_exomiser_file:$yml_exomiser_file,final_anno:$final_anno");
	return if($flag_dup == 1);


	#注释基因信息by snpeff
	&cmd_print($print_flag, "#&J	$task_id	$pre_task\n");
	&cmd_print($print_flag, "$bin_mkdir $data_dir");
	&cmd_print($print_flag, "$bin_annovar/table_annovar.pl -out $annovar_anno_out -buildver hg19 -protocol clinvar_20190305,cosmic70,civic_20200317,civicRegion4SmallMut_20200317,gwasCatalog,avsnp150,gnomad211_genome,gnomad211_exome,dbnsfp35c -operation f,f,f,r,r,f,f,f,f $final_vcf $annovar_humandb --remove --otherinfo --nastring . --vcfinput -intronhgvs 100000");

	return;
	&cmd_print("$bin_create_exomiser_exome liumm_vcf_file_path:$final_vcf liumm_output_prefix:$exomiser_anno_prefix >$yml_exomiser_file");
	#&cmd_print("$bin_snpeff -s $snpeff_anno_stat $final_vcf >$snpeff_anno_vcf");

	#my $bin_exomiser="java -Xms2g -Xmx4g -jar $bin_dir/../../software/exomiser-cli-12.1.0/exomiser-cli-12.1.0.jar";
	my $bin_create_exomiser_exome="cat $bin_dir/exomiser/exome.yml|perl $bin_dir/exomiser/create_yml.pl";

	&cmd_print("$bin_snpeff -s $snpeff_anno_stat $final_vcf >$snpeff_anno_vcf");
	&cmd_print("cat $snpeff_anno_vcf|$bin_fix_snpeff|$bin_select_snpeff|$bin_anno_gene -clinvar_all $database_clinvar_stat_all -clinvar_mul $database_clinvar_stat_mul -hpo_g2p $database_g2p_hpo  -hpo_p2g $database_p2g_hpo >$snpeff_anno_vcf.fix.select.annogene.vcf");
	#注释各类数据库信息by annovar
	#&cmd_print("$bin_annovar/table_annovar.pl -out $annovar_anno_out -buildver hg19 -protocol clinvar_20180603,cosmic70,civic_20181226,civicRegion_20181226,avsnp150,gnomad_genome,gnomad_exome,1000g2015aug_all,1000g2015aug_eas,dbnsfp33a -operation f,f,f,r,f,f,f,f,f,f $snpeff_anno_vcf.fix.select.annogene.vcf $annovar_humandb --remove --otherinfo --nastring . --vcfinput");
	&cmd_print("$bin_annovar/table_annovar.pl -out $annovar_anno_out -buildver hg19 -protocol clinvar_20190305,cosmic70,civic_20200317,civicRegion4SmallMut_20200317,gwasCatalog,avsnp150,gnomad211_genome,gnomad211_exome,dbnsfp35c -operation f,f,f,r,r,f,f,f,f $snpeff_anno_vcf.fix.select.annogene.vcf $annovar_humandb --remove --otherinfo --nastring . --vcfinput");
	#&cmd_print("$bin_annovar/table_annovar.pl -out $somatic_annovarOut -buildver hg19 -protocol refGeneWithVer,cosmic70,civic_20200317,civicRegion4SmallMut_20200317 -operation gx,f,f,r $somatic_avinput $annovar_humandb  --remove --otherinfo --nastring -");

}


sub vcf_stat{

	my $task_id = "AllStatVcf";

	my $pre_task = "";
	my $input_str = "";
	foreach my $sample_id( keys(%{$h_sample_info{"sample"}}) ){
		foreach my $metric_id('final_vcf'){
			if( defined($h_sample_info{"sample"}{$sample_id}{$metric_id}) ){
				my $anno_type = "germline";
				$input_str .= " $sample_id\_$anno_type:$anno_type:".$h_sample_info{"sample"}{$sample_id}{$metric_id}{"path"};
				$pre_task .= ",".$h_infile4preTask{$h_sample_info{"sample"}{$sample_id}{$metric_id}{"path"}};
			}
		}
	}

	if($pre_task =~ s/^,//){
		my $data_dir = "$out_dir/final_result";
		$input_str =~ s/^,//;
		&cmd_print($print_flag, "#&J	$task_id	$pre_task\n");
		&cmd_print($print_flag, "$bin_mkdir $data_dir");
		&cmd_print($print_flag, "$bin_smallVariant $bed_region $data_dir/AllStat.smallVariant_stat.detail.xls $data_dir/AllStat.smallVariant_stat.stat.xls $input_str");
	}
}




exit(0);












my $line;
open IN, $infile or die $!;
while( $line = <IN> ){
	next if( $line =~ /^\s*$/ );
	chomp $line;
	#mode	gender	sample1	s1_fq1	s1_fq2	sample2	s2_fq1	s2_fq2
	my @arr = split(/\t/,$line);

	#获取bed文件，pon文件等配置数据
	if( $arr[0] =~ /^#/ ){
		if( $arr[0] =~ /^#target_bed\s*$/i ){
			#获取bed文件
			$bed_region = $arr[1];
			#获取interval_list文件for GATK & Picard
			$gatk_interval_list = $arr[2];
			$h_infile4docker{$bed_region}++;
			$h_infile4docker{$gatk_interval_list}++;
		}
		elsif( $arr[0] =~ /^#pon4gatk\s*$/i ){
			#获取mutect2检测somatic变异需要用到的panel of normal数据
			$pon4gatk_vcf = $arr[1];
		}
		elsif( $arr[0] =~ /^#cnvkit_odds_limit_male\s*$/i ){
			#获取用于CNVkit的男性的log2 odds阈值
			$cnvkit_odds_limit_male = $arr[1];
		}
		elsif( $arr[0] =~ /^#cnvkit_odds_limit_female\s*$/i ){
			#获取用于CNVkit的女性的log2 odds阈值
			$cnvkit_odds_limit_female = $arr[1];
		}
		elsif( $arr[0] =~ /^#gatk_cnv_interval_list\s*$/i ){
			#获取用于GATK CNV的interval list文件
			$gatk_cnv_interval_list = $arr[1];
		}
		elsif( $arr[0] =~ /^#normal_match\s*$/i ){
			#获取已运行完使用过的normal样本作为对照
			#$recal_bam2	$sorted_bam2	$table_getpileupsummaries2
			shift@arr;
			my $normal_sample_id = shift@arr;
			#默认 Job_ID为NA
			push @arr, "NA";
			$h_normal{$normal_sample_id} = @arr;
		}
		next;
	}

	#获取测序数据信息
	my $t_modes = shift@arr;
	if( $t_modes =~ /sample_pe/i ){#1. PE测序原始数据处理
		if( $t_modes =~ /^\s*sample_pe\s*$/i ){#1. PE测序原始数据处理
			#获取样本信息
			my ( $gender, $sample,$hpo_ids, $fq1, $fq2 ) = @arr;
			if( defined($h_sample_info{"sample"}{$sample}) ){
				&err_print("Err: duplicate sample $sample");
				exit(0);
			}

			$h_sample_info{"sample"}{$sample}{"id"} = $sample;
			$h_sample_info{"sample"}{$sample}{"gender"} = $gender;
			$h_sample_info{"sample"}{$sample}{"raw_fq1"}{"path"} = $fq1;	$h_sample_info{"sample"}{$sample}{"raw_fq1"}{"task"} = "NA";
			$h_sample_info{"sample"}{$sample}{"raw_fq2"}{"path"} = $fq2;	$h_sample_info{"sample"}{$sample}{"raw_fq2"}{"task"} = "NA";
			$h_sample_info{"sample"}{$sample}{"clean_fq1"}{"path"} = $fq1;	$h_sample_info{"sample"}{$sample}{"clean_fq1"}{"task"} = "NA";
			$h_sample_info{"sample"}{$sample}{"clean_fq2"}{"path"} = $fq2;	$h_sample_info{"sample"}{$sample}{"clean_fq2"}{"task"} = "NA";
			$h_infile4docker{$fq1}++;
			$h_infile4docker{$fq2}++;
			#&fastp_qc($sample);
			&gatk_preprocess($sample);
			$print_flag = 1;
			&gatk_qc($sample);

			#获取HPO表型
			if( $hpo_ids ne "-" ){
				foreach my $hpo_info( split(/\|/,$hpo_ids) ){
					if( $hpo_info =~ /^\s*(HP:\d+)\s*-/i ){
						$h_sample_info{"sample"}{$sample}{"hpo"}{$1}++;
					}
					else{
						print STDERR "Err: illegal HPO format for $sample\n";
						exit(0);
					}
				}
			}

		}
		elsif( $t_modes =~ /^\s*create_pon4gatk\s*$/i ){#1. PE测序原始数据处理
			my ( $pon_id, $samples ) = @arr;
			&create_pon4gatk( $pon_id, $samples );
		}
		else{
			print STDERR "Err: illegal mode: $t_modes\n";
			exit(0);
		}
		next;
	}
	foreach my $mode( split(/,/,$t_modes) ){
		if( $mode =~ /^\s*pre_notrim\s*$/i ){
			#1. 比对，BQSR等
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&preprocess("notrim",$sample);
					}
				}
			}
		}
		elsif( $mode =~ /^\s*pre_trim\s*$/i ){
			#1. 比对，BQSR等
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&preprocess("trim",$sample);
					}
				}
			}
		}
		elsif( $mode =~ /^\s*QC\s*$/i ){
			#2. QC
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&quality_control( $sample );
					}
				}
			}
		}
		elsif( $mode =~ /^\s*HLA_typing\s*$/i ){
			#1. 比对，BQSR等
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&hla_typing( $sample );
					}
				}
			}
		}
		elsif( $mode =~ /^\s*factera\s*$/i ){
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&factera( $sample );
					}
				}
			}
		}
		elsif( $mode =~ /^\s*msisensor_pair\s*$/i ){#3.1 肿瘤Pair样本的MSI检测
			foreach my $samples( @arr ){
				my ( $sample_tumor, $sample_normal ) = split(/,/, $samples);
				if( $sample_tumor =~ s/^\s*(\S+)\s*$/$1/ && $sample_normal =~ s/^\s*(\S+)\s*$/$1/ ){
					&msisensor_pair( $sample_tumor, $sample_normal );
				}
			}
		}
		elsif( $mode =~ /^\s*msisensor_single\s*$/i ){#3.2 肿瘤Single样本的MSI检测
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&msisensor_single($sample);
					}
				}
			}
		}
		elsif( $mode =~ /^\s*cnvkit_pair\s*$/i ){#4. Pair样本的CNV检测
			foreach my $samples( @arr ){
				my ( $sample_tumor, $sample_normal ) = split(/,/, $samples);
				if( $sample_tumor =~ s/^\s*(\S+)\s*$/$1/ && $sample_normal =~ s/^\s*(\S+)\s*$/$1/ ){
					&cnvkit_pair( $sample_tumor, $sample_normal );
				}
			}
		}
		elsif( $mode =~ /^\s*germline_gatk4\s*$/i ){#5. Single样本的Germline变异检测
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&germline_gatk4( $sample );
					}
				}
			}
		}
		elsif( $mode =~ /^\s*tumor_pair_gatk4\s*$/i ){#总2：肿瘤pair样本分析
			#my ( $sample, $normal_sample ) = @arr;
			foreach my $samples( @arr ){
				my ( $sample_tumor, $sample_normal ) = split(/,/, $samples);
				if( $sample_tumor =~ s/^\s*(\S+)\s*$/$1/ && $sample_normal =~ s/^\s*(\S+)\s*$/$1/ ){
					&somatic_callVairant( "gatk4", $sample_tumor, $sample_normal);
				}
			}
		}
		elsif( $mode =~ /^\s*tumor_pair_vardict\s*$/i ){#总2：肿瘤pair样本分析
			#my ( $sample, $normal_sample ) = @arr;
			foreach my $samples( @arr ){
				my ( $sample_tumor, $sample_normal ) = split(/,/, $samples);
				if( $sample_tumor =~ s/^\s*(\S+)\s*$/$1/ && $sample_normal =~ s/^\s*(\S+)\s*$/$1/ ){
					&somatic_callVairant( "vardict", $sample_tumor, $sample_normal);
				}
			}
		}
		elsif( $mode =~ /^\s*tumor_only_gatk4\s*$/i ){#总2：肿瘤pair样本分析
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&somatic_callVairant( "gatk4", $sample, "NoNormal");
					}
				}
			}
		}
		elsif( $mode =~ /^\s*tumor_pon_gatk4\s*$/i ){#总2：肿瘤pair样本分析
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&somatic_callVairant( "gatk4", $sample, "PON");
					}
				}
			}
		}
		elsif( $mode =~ /^\s*tumor_only_vardict\s*$/i ){#总2：肿瘤pair样本分析
			foreach my $samples( @arr ){
				foreach my $sample( split(/,/,$samples) ){
					if($sample =~ s/^\s*(\S+)\s*$/$1/){
						&somatic_callVairant( "vardict", $sample, "NoNormal");
					}
				}
			}
		}
		else{
			&err_print("Err: unrecognized mode $mode");
			exit(0);
		}
	}
}
close IN;

#综合统计
&anno_all_variant();#注释所有变异检测结果
&sex_stat();
&msi_stat();
&factera_stat();
&QC_stat();
&HLA_stat();


sub task_adminbak{
	
	#初始化
	my ($task_id, $sample_id, $in_files, $out_files)= @_;
	my %pre_ids = ("NA"=>1);
	my @arr_return = (0);


	#判断该任务是否已运行过
	if( defined( $h_cmd_log{$task_id} ) ){
		$arr_return[0] = 1;
		return @arr_return;
	}
	$h_cmd_log{$task_id}++;

	#输入文件处理
	if( $in_files ne "-" ){
		foreach my $file( split(/,/,$in_files) ){
			push @arr_return,$h_sample_info{"sample"}{$sample_id}{$file}{"path"};
			$pre_ids{ $h_sample_info{"sample"}{$sample_id}{$file}{"task"} } ++;
		}
	}
	my $str_pre_id = join(",",sort {$a cmp $b} keys(%pre_ids));
	$str_pre_id =~ s/,NA,//;	$str_pre_id =~ s/^NA,//;	$str_pre_id =~ s/,NA$//;
	push @arr_return, $str_pre_id;

	#输入文件处理
	if( $out_files ne "-" ){
		foreach my $name2file( split(/,/,$out_files) ){
			my ($name,$file) = split(/:/,$name2file);
			$h_sample_info{"sample"}{$sample_id}{$name}{"path"} = $file;
			$h_sample_info{"sample"}{$sample_id}{$name}{"task"} = $task_id;
		}
	}


	return @arr_return;
}


#Required by GATK CNV pipeline
sub gatk_cnvIntervals{
	#my ( $data_dir ) = @_;
	#任务ID
	my $task_id = "gatk_interval";
	my $data_dir="$out_dir/gatk_interval";

	#&cmd_print("$bin_gatk PreprocessIntervals -R $ref_fa -L $bed_region --bin-length 0 --interval-merging-rule OVERLAPPING_ONLY -O out_cnv_interval_list");
	#&cmd_print("$bin_gatk PreprocessIntervals -R $ref_fa -L $bed_region --bin-length 0 --interval-merging-rule OVERLAPPING_ONLY -O $cnv_interval_list");
}


#1. preprocess
sub fastq_qcbak{
	
	#任务ID
	my ($sample_id) = @_;
	my $task_id = "$sample_id\_fastq_QC";
	#初始化文件名
	my $data_dir="$out_dir/fastp/$sample_id";
	my $trimFq_prefix="$data_dir/$sample_id.trim";
	my $out_trimFq1 = "$trimFq_prefix.R1.fq.gz";
	my $out_trimFq2 = "$trimFq_prefix.R2.fq.gz";

	my ($flag_dup, $raw_fq1, $raw_fq2, $pre_task)  = &task_admin($task_id, $sample_id, "raw_fq1,raw_fq2", "clean_fq1:$out_trimFq1,clean_fq2:$out_trimFq2");
	return if($flag_dup == 1);
	
	#运行脚本
	&cmd_print("#&J	$task_id	$pre_task");
	&cmd_print("$bin_mkdir $data_dir ");
	&cmd_print("$bin_fastp -i $raw_fq1 -I $raw_fq2 -o $trimFq_prefix.R1.fq.gz -O $trimFq_prefix.R2.fq.gz  -h $trimFq_prefix.html -j $trimFq_prefix.json");

}






#1. preprocess
sub preprocess{

	#args
	my ($trim_flag,$sample_id) = @_;
	my $umi_flag = "no";
	if( defined( $h_cmd_log{"preprocess"}{$sample_id} ) ){
		return;
	}
	else{
		$h_cmd_log{"preprocess"}{$sample_id} = 1;
	}

	#get input
	my $raw_fq1 = $h_sample_info{"sample"}{$sample_id}{"raw_fq1"};
	my $raw_fq2 = $h_sample_info{"sample"}{$sample_id}{"raw_fq2"};
	
	#set names
	my $data_dir="$out_dir/preprocess/$sample_id";     `mkdir -p $data_dir`;
	my $trimFq_prefix="$data_dir/$sample_id.trim";
	my $bwa_prefix="$data_dir/$sample_id";
	my $raw_sam = "$bwa_prefix\.sam";
	my $markDup_bam = "$bwa_prefix\.markDup.bam";
	my $markDup_metric = "$bwa_prefix\.markDup.metrics.txt";
	my $sorted_bam="$bwa_prefix\.markDup.sorted.bam";	$h_sample_info{"sample"}{$sample_id}{"sorted_bam"} = $sorted_bam;
	my $recal_table="$bwa_prefix\.markDup.sorted.recal_data.table";#$h_sample_info{"sample"}{$sample_id}{"recal_table"} = $recal_table;
	my $recal_bam="$bwa_prefix\.markDup.sorted.recal.bam";	$h_sample_info{"sample"}{$sample_id}{"recal_bam"} = $recal_bam;
	my $table_getpileupsummaries="$recal_bam\.getpileupsummaries.table";
	$h_sample_info{"sample"}{$sample_id}{"table_getpileupsummaries"} = $table_getpileupsummaries;

	#preprocess following GATK best practice
	{
		&cmd_print("#&J	$sample_id\_aln	NA");
		my $in_fq1 = $raw_fq1;
		my $in_fq2 = $raw_fq2;
		#判读是否过滤
		if( $trim_flag =~ /^trim$/i ){
			&cmd_print("$bin_trimmomatic PE -phred33 -threads 8 -summary $trimFq_prefix\.summary.txt $raw_fq1 $raw_fq2 $trimFq_prefix\_PE_R1.fq.gz $trimFq_prefix\_SE_R1.fq.gz $trimFq_prefix\_PE_R2.fq.gz $trimFq_prefix\_SE_R2.fq.gz MINLEN:75");

			$in_fq1 = "$trimFq_prefix\_R1.fq.gz";
			$in_fq2 = "$trimFq_prefix\_R2.fq.gz";
		}
		#判读是否进行UMI处理
		if( $umi_flag !~ /^noa$/i ){
			my $filterumi_fq1 = "$data_dir/$sample_id\.filter_umi_R1.fq";
			my $filterumi_fq2 = "$data_dir/$sample_id\.filter_umi_R2.fq";
			my $filterumi_stat = "$data_dir/$sample_id\.filter_umi.stat.xls";
			my $unmap_mark_umi_bam = "$data_dir/$sample_id\.unmap.mark_umi.bam";
			my $unmap_mark_umi_fq = "$data_dir/$sample_id\.unmap.mark_umi.fq";
			my $map_mark_umi_bam = "$data_dir/$sample_id\.map.mark_umi.bam";
			my $all_mark_umi_bam = "$data_dir/$sample_id\.all.mark_umi.bam";
			my $group_mark_umi_bam = "$data_dir/$sample_id\.group.mark_umi.bam";
			my $concensus_mark_umi_bam = "$data_dir/$sample_id\.concensus.mark_umi.bam";
			my $concensus_mark_umi_fq = "$data_dir/$sample_id\.concensus.mark_umi.fq";
			my $concensus_map_bam = "$data_dir/$sample_id\.concensus.map.bam";
			my $concensus_all_bam = "$data_dir/$sample_id\.concensus.all.bam";
			my $concensus_filter_bam = "$data_dir/$sample_id\.concensus.filter.bam";
			my $concensus_filter_sort_bam = "$data_dir/$sample_id\.concensus.filter.sorted.bam";
			my $concensus_clip_bam = "$data_dir/$sample_id\.concensus.clip.bam";
			&cmd_print("$bin_filter_umi_cfbestV2 $in_fq1 $in_fq2 $filterumi_fq1 $filterumi_fq2 $filterumi_stat");
			&cmd_print("$bin_fgbio FastqToBam --sort true --sample $sample_id --library lib_$sample_id --input $filterumi_fq1 $filterumi_fq2 --read-structures 7M1S+T 7M1S+T --output $unmap_mark_umi_bam");
			&cmd_print("$bin_picard SamToFastq I=$unmap_mark_umi_bam F=$unmap_mark_umi_fq INTERLEAVE=true");
			&cmd_print("$bin_bwa mem -t 16 -M -p $ref_fa $unmap_mark_umi_fq >$map_mark_umi_bam");
			&cmd_print("$bin_picard MergeBamAlignment ALIGNED=$map_mark_umi_bam UNMAPPED=$unmap_mark_umi_bam R=$ref_fa O=$all_mark_umi_bam");
			&cmd_print("$bin_fgbio GroupReadsByUmi --input $all_mark_umi_bam --output $group_mark_umi_bam --strategy=paired --min-map-q=20 --edits=1");
			&cmd_print("$bin_fgbio CallMolecularConsensusReads --input $group_mark_umi_bam --output $concensus_mark_umi_bam --min-reads=1 --min-input-base-quality=20");
			&cmd_print("$bin_fgbio FilterConsensusReads --input $concensus_mark_umi_bam --output $concensus_filter_bam --ref $ref_fa --min-reads=1 --max-read-error-rate=0.05 --min-base-quality=20");
			&cmd_print("$bin_picard SortSam I=$concensus_filter_bam SORT_ORDER=queryname O=$concensus_filter_sort_bam");
			&cmd_print("$bin_picard SamToFastq I=$concensus_filter_sort_bam F=$concensus_mark_umi_fq INTERLEAVE=true");
			&cmd_print("$bin_bwa mem -t 16 -M -p $ref_fa $concensus_mark_umi_fq >$concensus_map_bam");
			&cmd_print("$bin_picard MergeBamAlignment ALIGNED=$concensus_map_bam UNMAPPED=$concensus_filter_sort_bam R=$ref_fa O=$concensus_all_bam");
			#&cmd_print("$bin_fgbio FilterConsensusReads --input $concensus_all_bam --output $concensus_filter_bam --ref $ref_fa --min-reads=1 --max-read-error-rate=0.05 --min-base-quality=20");
			#&cmd_print("$bin_fgbio ClipBam --input $concensus_all_bam --output $concensus_clip_bam --ref $ref_fa --clip-overlapping-reads=true");
			$markDup_bam = $concensus_all_bam;
			$h_sample_info{"sample"}{$sample_id."-nofilter"}{"sorted_bam"} = $all_mark_umi_bam;
			&quality_control( $sample_id."-nofilter" );
		}
		else{
			&cmd_print("$bin_bwa mem -t 16 -M $ref_fa $in_fq1 $in_fq2 >$raw_sam");
			&cmd_print("$bin_picard MarkDuplicates I=$raw_sam O=$markDup_bam M=$markDup_metric VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname");
			$h_sample_info{"sample"}{$sample_id}{"QC"}{$markDup_metric}++;
		}
		&cmd_print("$bin_picard AddOrReplaceReadGroups I=$markDup_bam O=$sorted_bam SORT_ORDER=coordinate RGID=$sample_id  RGLB=lb_$sample_id RGPL=illumina RGPU=pu_$sample_id RGSM=$sample_id TMP_DIR=$temp_dir");
		&cmd_print("$bin_samtools index $sorted_bam");
		&cmd_print("$bin_gatk BaseRecalibrator -R $ref_fa -I $sorted_bam -O $recal_table --known-sites $gatk_dbsnp --known-sites $gatk_Mills_and_1000G_gold_standard_indels --use-original-qualities");
		&cmd_print("$bin_gatk ApplyBQSR -R $ref_fa -I $sorted_bam -bqsr $recal_table -O $recal_bam  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30  --add-output-sam-program-record --create-output-bam-md5  --use-original-qualities");
		&cmd_print("$bin_gatk GetPileupSummaries -R $ref_fa -V $gatk_GetPileupSummaries --interval-set-rule INTERSECTION -L $gatk_interval_list -L $gatk_GetPileupSummaries -I $recal_bam -O $table_getpileupsummaries");

		#不用GATK CNV流程
		#GATK CNV required
		#&cmd_print("#&J	$sample_id\_CollectReadCounts	$sample_id\_aln");
		#&cmd_print("$bin_gatk CollectReadCounts -R $ref_fa -I $recal_bam -L $gatk_cnv_interval_list --interval-merging-rule OVERLAPPING_ONLY -O $read_count --format TSV");
		#&cmd_print("$bin_gatk CollectAllelicCounts -R $ref_fa -I $recal_bam -L $bed_region -O $allelic_counts");
	}
}

#HLA typing
sub hla_typing{

	my ($sample_id) = @_;;
	if( defined( $h_cmd_log{"HLA_Typing"}{$sample_id} ) ){
		return;
	}
	else{
		$h_cmd_log{"HLA_Typing"}{$sample_id} = 1;
	}


	#get input
	my $sorted_bam = $h_sample_info{"sample"}{$sample_id}{"sorted_bam"};
	my $dir_bam=$out_dir;
	my $file_bam=$sorted_bam;
	if($file_bam =~ s/^$out_dir\/+([^\/])/$1/){
	}
	else{
		print STDERR "illegal bam path $sorted_bam\n";
		exit(0);
	}
	#set names
	my $data_dir="$out_dir/HLA/$sample_id";	`mkdir -p $data_dir`;
	my $qc_prefix="$data_dir/$sample_id.qc";

	#QC by picard & bedtools
	my $last_job = "$sample_id\_aln";
	&cmd_print("#&J	$sample_id\_HLA	$last_job");
	&cmd_print("NAME=\"polysolver_$sample_id\";	DIR1=\"$dir_bam\";	DIR2=\"/home/docker/polysolver_$sample_id\";	OUT_DIR=\"HLA/$sample_id\";	BAM=\"$file_bam\"");
	&cmd_print("$bin_polysolver");
	$h_sample_info{"sample"}{$sample_id}{"HLA"} = $data_dir."/genome_results.txt";
	#my $bin_polysolver="/usr/bin/docker run -P --name \$NAME -v \$DIR1:\$DIR2 sachet/polysolver:v4 bash /home/polysolver/scripts/shell_call_hla_type \$DIR2/\$BAM Unknown 1 hg19 STDFQ 0 \$DIR2/\$OUT_DIR; /usr/bin/docker rm $NAME";
}


#QC
sub quality_control{

	my ($sample_id) = @_;;
	if( defined( $h_cmd_log{"quality_control"}{$sample_id} ) ){
		return;
	}
	else{
		$h_cmd_log{"quality_control"}{$sample_id} = 1;
	}


	#get input
	my $sorted_bam = $h_sample_info{"sample"}{$sample_id}{"sorted_bam"};
	#set names
	my $data_dir="$out_dir/QC/$sample_id";	`mkdir -p $data_dir`;
	my $qc_prefix="$data_dir/$sample_id.picardQC";

	#QC by picard & bedtools
	my $last_job = "$sample_id\_aln";
	&cmd_print("#&J	$sample_id\_QC	$last_job");
	my $t_interval1 = "";
	my $t_interval2 = "";
	if( defined($gatk_interval_list) && $gatk_interval_list =~ /\S/ ){
		$t_interval1 = "INTERVALS=$gatk_interval_list";
		$t_interval2 = "BAIT_INTERVALS=$gatk_interval_list  TARGET_INTERVALS=$gatk_interval_list";
	}
	#&cmd_print("$bin_qualimap bamqc -bam $sorted_bam -gff $bed_region -outdir $data_dir -outformat pdf");
	#$h_sample_info{"sample"}{$sample_id}{"QC"} = $data_dir."/genome_results.txt";
	$h_sample_info{"sample"}{$sample_id}{"QC"}{"$qc_prefix*metrics*"}++;
	my $program = "PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectGcBiasMetrics PROGRAM=CollectSequencingArtifactMetrics PROGRAM=CollectQualityYieldMetrics";;
	&cmd_print("$bin_picard CollectMultipleMetrics $program I=$sorted_bam R=$ref_fa O=$qc_prefix.multiple_metrics $t_interval1");
	&cmd_print("$bin_picard CollectHsMetrics I=$sorted_bam R=$ref_fa O=$qc_prefix.hs_metrics.txt $t_interval2 PER_TARGET_COVERAGE=$qc_prefix.hs_metrics.per_target.txt COVERAGE_CAP=50000");
	#&cmd_print("sed '1,11d' $qc_prefix.hs_metrics.txt|awk '\$2>0' > $qc_prefix.hs_metrics.plot_data.txt");
	#&opt_print("$bin_picard CollectHsMetrics I=$sorted_bam R=$ref_fa O=$qc_prefix.hs_metrics.mq0_bq0.txt $t_interval2 PER_TARGET_COVERAGE=$qc_prefix.hs_metrics.per_target.mq0_bq0.txt COVERAGE_CAP=50000 MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0");
	#&opt_print("$bin_picard CollectWgsMetrics I=$sorted_bam R=$ref_fa O=$qc_prefix.wgs_metrics.txt $t_interval1 COVERAGE_CAP=50000 MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 COUNT_UNPAIRED=true");
	#&opt_print("$bin_picard CollectGcBiasMetrics I=$sorted_bam R=$ref_fa O=$qc_prefix.gc_bias_metrics.txt S=$qc_prefix.gc_summary_metrics.txt CHART=$qc_prefix.gc_bias_metrics.pdf ALSO_IGNORE_DUPLICATES=true");
	#bedtools coverage to calculate coverage
	#&cmd_print("$bin_samtools view -bF 1284 $sorted_bam|$bin_bedtools coverage -hist -a $bed_region -b stdin >$qc_prefix.bed_cov.txt");
	#&cmd_print("grep '^all' $qc_prefix.bed_cov.txt|awk -F \"\\t\" '{print \$2\"\\t\"\$3}' >$qc_prefix.bed_cov.plot_data.txt");
}


#fusion detection with factera
sub factera{

	my ($sample_id) = @_;
	if( defined( $h_cmd_log{"factera"}{$sample_id} ) ){
		return;
	}
	else{
		$h_cmd_log{"factera"}{$sample_id} = 1;
	}

	#get input
	my $sorted_bam = $h_sample_info{"sample"}{$sample_id}{"sorted_bam"};
	my $sorted_bam_prefix = basename($sorted_bam);
	$sorted_bam_prefix =~ s/.bam$//;

	#set names
	my $data_dir="$out_dir/factera/$sample_id";	`mkdir -p $data_dir`;
	my $factera_out_dir = "$data_dir/";
	my $factera_out = "$factera_out_dir/$sorted_bam_prefix\.factera.fusions.txt";	$h_sample_info{"sample"}{$sample_id}{"factera_out"} = $factera_out;

	#factera SV(e.g. fusion)
	my $last_job = "$sample_id\_aln";
	&cmd_print("#&J	$sample_id\_factera	$last_job");
	my $t_bed = "";
	$t_bed  = $bed_region if( defined($bed_region) && $bed_region =~ /\S/ );
	&cmd_print("$bin_factera -C -o $factera_out_dir $sorted_bam $factera_bed $ref_fa_2bit $t_bed");

}

sub msisensor_single{

	my  ( $sample_id ) = @_;
	if( defined( $h_cmd_log{"msisensor_single"}{$sample_id} ) ){
		return;
	}
	else{
		$h_cmd_log{"msisensor_single"}{$sample_id} = 1;
	}

	#get input
	my $sorted_bam = $h_sample_info{"sample"}{$sample_id}{"sorted_bam"};
	#set names
	my $data_dir="$out_dir/msisensor/$sample_id";	`mkdir -p $data_dir`;
	my $msisensor_out = "$data_dir/$sample_id\.single.msi";

	#MSI by msisensor
	my $last_job = "$sample_id\_aln";
	&cmd_print("#&J	$sample_id\_msisensor	$last_job");
	&cmd_print("$bin_msisensor msi -d $ref_msi_microsatellite -t $sorted_bam -e $bed_region -o $msisensor_out");
	$h_sample_info{"sample"}{"$sample_id"}{"msi"} = $msisensor_out;

}

sub msisensor_pair{
	my  ($sample_id, $normal_sample_id ) = @_;
	if( defined( $h_cmd_log{"msisensor_pair"}{"$sample_id\_VS_$normal_sample_id"} ) ){
		return;
	}
	else{
		$h_cmd_log{"msisensor_pair"}{"$sample_id\_VS_$normal_sample_id"} = 1;
	}

	#get input
	my $sorted_bam = $h_sample_info{"sample"}{$sample_id}{"sorted_bam"};
	my $normal_sorted_bam = $h_sample_info{"sample"}{$normal_sample_id}{"sorted_bam"};
	#set names
	my $data_dir = "$out_dir/msisensor/$sample_id\_VS_$normal_sample_id";	`mkdir -p $data_dir`;
	my $msisensor_out = "$data_dir/$sample_id\_VS_$normal_sample_id\.msi";

	#MSI by msisensor
	my $last_job = "$sample_id\_aln,$normal_sample_id\_aln";
	&cmd_print("#&J	$sample_id\_VS_$normal_sample_id\_msisensor	$last_job");
	&cmd_print("$bin_msisensor msi -d $ref_msi_microsatellite -n $normal_sorted_bam -t $sorted_bam -e $bed_region -o $msisensor_out");
	$h_sample_info{"sample"}{"$sample_id\_VS_$normal_sample_id"}{"msi"} = $msisensor_out;

}

#Call CNV by pair
sub cnvkit_pair{
	my  ($sample_id, $normal_sample_id ) = @_;
	if( defined( $h_cmd_log{"cnvkit_pair"}{"$sample_id\_VS_$normal_sample_id"} ) ){
		return;
	}
	else{
		$h_cmd_log{"cnvkit_pair"}{"$sample_id\_VS_$normal_sample_id"} = 1;
	}

	#get input
	my $sorted_bam = $h_sample_info{"sample"}{$sample_id}{"sorted_bam"};
	my $sorted_bam_prefix = basename($sorted_bam);
	$sorted_bam_prefix =~ s/.bam$//;
	my $normal_sorted_bam = $h_sample_info{"sample"}{$normal_sample_id}{"sorted_bam"};
	#set names
	my $data_dir="$out_dir/cnvkit/$sample_id\_VS_$normal_sample_id";	`mkdir -p $data_dir`;
	my $cnv_out_dir = "$data_dir";
	my $cnvkit_outref="$data_dir/$sample_id\_VS_$normal_sample_id\.cnvkit.reference.cnn";

	#cnvkit
	my $last_job = "$sample_id\_aln,$normal_sample_id\_aln";
	&cmd_print("#&J	$sample_id\_VS_$normal_sample_id\_cnvkit	$last_job");
	&cmd_print("$bin_cnvkits batch $sorted_bam --normal $normal_sorted_bam --targets $bed_region --annotate $ref_flat --fasta $ref_fa --access $ref_cnvkits_access --output-reference $cnvkit_outref --output-dir $cnv_out_dir --diagram --scatter");
	
	my $cnv_odds = $cnvkit_odds_limit_common;
	#获取对应性别的CNV过滤标准
	my $gender = $h_sample_info{"sample"}{$sample_id}{"gender"};
	if( $gender =~ /^f/i ){
		$cnv_odds = $cnvkit_odds_limit_female;
	}
	elsif( $gender =~ /^m/i ){
		$cnv_odds = $cnvkit_odds_limit_male;
	}
	else{
		&err_print("Warn: illegal gender $gender for $sample_id. Default cnvkit_odds_limit: $cnv_odds");
	}
	&cmd_print("$bin_filterCNVkit -odds $cnv_odds < $cnv_out_dir/$sorted_bam_prefix\.cns >$cnv_out_dir/$sorted_bam_prefix\.filtered.cns");
             #perl $0 -odds 0.3,-0.4,0.3,-0.4,0.3,-0.4 <$cnv_out_dir/$sorted_bam_prefix\.cns  >B1701.markDup.sorted.filtered.cns
	#&cmd_print("$bin_filterCNVkit -odds $cnv_odds -limit batch $bam_file $nomal_cnvkit_cmd --targets $bed_region --annotate $ref_flat --fasta $ref_fa --access $ref_cnvkits_access --output-reference $cnvkit_outref --output-dir $data_dir --diagram --scatter");

}




#检测任务是否重复
sub check_dup_job{
	my  ($job_name, $sample_name) = @_;
	if( defined( $h_cmd_log{$job_name}{$sample_name} ) ){
		print STDERR "Warning:duplicated job:$job_name\t$sample_name\n";
		return 1;
	}
	else{
		$h_cmd_log{$job_name}{$sample_name} = 1;
		return 0;
	}
}

sub create_pon4gatk{
	my ($pon_id,$sample_ids) = @_;

	my %h_temp_samples = ();
	my $data_dir="$out_dir/pon4gatk";	`mkdir -p $data_dir`;
	my $pon_db = "$data_dir/pon_db_$pon_id";
	my $pon_vcf = "$data_dir/$pon_id\.pon4gatk.vcf";
	my $last_job = "";
	my $mutect_cmd = "";
	foreach my $sample_id( split(/,/,$sample_ids) ){
		my $recal_bam  = $h_sample_info{"sample"}{$sample_id}{"recal_bam"};
		my $somatic_raw_vcf = "$data_dir/$sample_id\.normalraw.vcf";
		if( defined($h_temp_samples{$somatic_raw_vcf}) ){
			print STDERR "Warning:duplicated sample for pon4gatk:$sample_id\n";
			next;
		}
		$last_job .= ",$sample_id\_aln";
		$mutect_cmd .= "$bin_gatk Mutect2 -R $ref_fa -I $recal_bam --max-mnp-distance 0 -O $somatic_raw_vcf\n";
		$h_temp_samples{$somatic_raw_vcf}++;
	}
	$last_job =~ s/^,//;
	my $input_vcfs = join(" -V ",sort {$a cmp $b} keys(%h_temp_samples));
	&cmd_print("#&J	$pon_id\_CreatePon4GATK	$last_job");
	&cmd_print("$mutect_cmd");
	&cmd_print("$bin_gatk GenomicsDBImport -R $ref_fa -L $gatk_interval_list --genomicsdb-workspace-path $pon_db -V $input_vcfs");
	&cmd_print("$bin_gatk CreateSomaticPanelOfNormals -R $ref_fa --germline-resource $mutect_gnomad -V gendb://$pon_db -O $pon_vcf");
}

sub somatic_callVairant{

	my  ($call_type, $sample_id, $normal_sample_id ) = @_;

	#get input
	my $sorted_bam = $h_sample_info{"sample"}{$sample_id}{"sorted_bam"};
	my $recal_bam  = $h_sample_info{"sample"}{$sample_id}{"recal_bam"};
	#set names
	#my $data_dir="$out_dir/somatic_call/$sample_id";	`mkdir -p $data_dir`;

	#call somatic variant for single tumor with vardict
	if( $call_type =~ /^vardict$/i && $normal_sample_id =~ /NoNormal/i ){
		#不允许重复任务
		return if( &check_dup_job("somatic_callVairant_vardict","$sample_id\_VS_$normal_sample_id") );
		my $data_dir="$out_dir/somatic_call/vardict/$sample_id\_VS_$normal_sample_id";	`mkdir -p $data_dir`;
		#vardict output
		my $vardict_vcf = "$data_dir/vardict.$sample_id\_VS_$normal_sample_id.vcf";
		my $last_job = "$sample_id\_aln";
		&cmd_print("#&J	$sample_id\_CallVariantSomaticVarDict	$last_job");
		&cmd_print("$bin_vardict -G $ref_fa -f 0.001 -N $sample_id -b $sorted_bam -c 1 -S 2 -E 3 -g 4 $bed_region|$bin_vardict_teststrandbias|$bin_vardict_var2vcf_valid  -N $sample_id -E -f 0.001 >$vardict_vcf");

		#记录最终vcf文件
		$h_sample_info{"call_variant"}{"SomaticVarDict"}{"$sample_id\_VS_$normal_sample_id"}{"final_vcf"} = $vardict_vcf;
		#记录该检测所有hpo的表型
		if( defined($h_sample_info{"sample"}{$sample_id}{"hpo"}) ){
			foreach my $hpo_id( keys(%{$h_sample_info{"sample"}{$sample_id}{"hpo"}})){
				$h_sample_info{"call_variant"}{"SomaticVarDict"}{"$sample_id\_VS_$normal_sample_id"}{"hpo"}{$hpo_id}++;
			}
		}
	}

	if( $call_type =~ /^gatk4$/i ){

		#获取pon信息
		my $cmd_pon = "";
		if( $normal_sample_id =~ /^PON$/i ){
			if( $pon4gatk_vcf =~ /\S/ ){
				$cmd_pon = "-pon $pon4gatk_vcf";
			}
			else{
				print STDERR "Err: No pon specified for run_mode tumor_only_gatk4 with pon\n";
				exit;
			}
		}
		return if( &check_dup_job("somatic_callVairant_gatk4", "$sample_id\_VS_$normal_sample_id") );
		#gatk4 output
		my $table_getpileupsummaries = $h_sample_info{"sample"}{$sample_id}{"table_getpileupsummaries"};
		my $data_dir="$out_dir/somatic_call/gatk4/$sample_id\_VS_$normal_sample_id";	`mkdir -p $data_dir`;
		my $somatic_prefix = "$data_dir/mutect.$sample_id\_VS_$normal_sample_id";
		my $somatic_raw_vcf="$somatic_prefix\.raw.vcf";
		my $somatic_f1r2="$somatic_prefix\.f1r2.tar.gz";
		my $somatic_rom="$somatic_prefix\.read-orientation-model.tar.gz";
		my $somatic_bam="$somatic_prefix\.bam";
		my $table_calculatecontamination="$somatic_prefix\.calculatecontamination.table";
		my $table_segmentation="$somatic_prefix\.segmentation.table";
		my $somatic_oncefiltered_vcf="$somatic_prefix\.oncefiltered.vcf";

		#call somatic SNV & InDel by Mutect2
		if( $normal_sample_id =~ /NoNormal/i || $normal_sample_id =~ /^PON$/i ){#call somatic SNV & InDel by Mutect2 with tumor-only model
			my $last_job = "$sample_id\_aln";
			&cmd_print("#&J	$sample_id\_CallVariantSomaticGATK4	$last_job");
			&cmd_print("$bin_gatk Mutect2 -R $ref_fa -I $recal_bam -tumor $sample_id $cmd_pon --germline-resource $mutect_gnomad --af-of-alleles-not-in-resource 0.0000025 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -L $gatk_interval_list -O $somatic_raw_vcf --f1r2-tar-gz $somatic_f1r2 -bamout $somatic_bam");
			&cmd_print("$bin_gatk CalculateContamination -I $table_getpileupsummaries -O $table_calculatecontamination --tumor-segmentation $table_segmentation");
			#记录该检测所有hpo的表型
			if( defined($h_sample_info{"sample"}{$sample_id}{"hpo"}) ){
				foreach my $hpo_id( keys(%{$h_sample_info{"sample"}{$sample_id}{"hpo"}})){
					$h_sample_info{"call_variant"}{"SomaticGATK4"}{"$sample_id\_VS_$normal_sample_id"}{"hpo"}{$hpo_id}++;
				}
			}
		}
		else{#call somatic SNV & InDel by Mutect2 tumor-pair model
			my $normal_recal_bam = $h_sample_info{"sample"}{$normal_sample_id}{"recal_bam"};
			my $normal_table_getpileupsummaries = $h_sample_info{"sample"}{$normal_sample_id}{"table_getpileupsummaries"};
			my $last_job = "$sample_id\_aln,$normal_sample_id\_aln";
			&cmd_print("#&J	$sample_id\_CallVariantSomaticGATK4	$last_job");
			&cmd_print("$bin_gatk Mutect2 -R $ref_fa -I $recal_bam -tumor $sample_id -I $normal_recal_bam -normal $normal_sample_id $cmd_pon --germline-resource $mutect_gnomad --af-of-alleles-not-in-resource 0.0000025 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -L $gatk_interval_list -O $somatic_raw_vcf --f1r2-tar-gz $somatic_f1r2 -bamout $somatic_bam");
			&cmd_print("$bin_gatk CalculateContamination -I $table_getpileupsummaries -O $table_calculatecontamination --tumor-segmentation $table_segmentation -matched $normal_table_getpileupsummaries");
		}
		&cmd_print("$bin_gatk LearnReadOrientationModel -I $somatic_f1r2 -O $somatic_rom");
		&cmd_print("$bin_gatk FilterMutectCalls -R $ref_fa -V $somatic_raw_vcf --contamination-table $table_calculatecontamination -O $somatic_oncefiltered_vcf --tumor-segmentation $table_segmentation --ob-priors $somatic_rom");

		#记录最终vcf文件
		$h_sample_info{"call_variant"}{"SomaticGATK4"}{ "$sample_id\_VS_$normal_sample_id" }{"final_vcf"} = $somatic_oncefiltered_vcf;
		#记录该检测所有hpo的表型
		if( defined($h_sample_info{"sample"}{$sample_id}{"hpo"}) ){
			foreach my $hpo_id( keys(%{$h_sample_info{"sample"}{$sample_id}{"hpo"}})){
				$h_sample_info{"call_variant"}{"SomaticGATK4"}{"$sample_id\_VS_$normal_sample_id"}{"hpo"}{$hpo_id}++;
			}
		}
		if( defined($h_sample_info{"sample"}{$normal_sample_id}{"hpo"}) ){
			foreach my $hpo_id( keys(%{$h_sample_info{"sample"}{$normal_sample_id}{"hpo"}})){
				$h_sample_info{"call_variant"}{"SomaticGATK4"}{"$sample_id\_VS_$normal_sample_id"}{"hpo"}{$hpo_id}++;
			}
		}
	}
}

sub somatic_pair_varscan{

	my  ($sample_id, $normal_sample_id ) = @_;
	if( defined( $h_cmd_log{"somatic_varscan"}{"$sample_id\_VS_$normal_sample_id"} ) ){
		return;
	}
	else{
		$h_cmd_log{"somatic_varscan"}{"$sample_id\_VS_$normal_sample_id"} = 1;
	}

	#get input
	my $sorted_bam = $h_sample_info{"sample"}{$sample_id}{"sorted_bam"};
	my $normal_sorted_bam = $h_sample_info{"sample"}{$normal_sample_id}{"sorted_bam"};
	#set names
	my $data_dir="$out_dir/varscan";	`mkdir -p $data_dir`;
	my $sorted_pileup="$sorted_bam\.mpileup";
	my $normal_sorted_pileup="$normal_sorted_bam\.mpileup";
	my $varscan_somatic_raw_vcf="$data_dir/varscan.$sample_id\_VS_$normal_sample_id\.raw.vcf";

	#mpileup required by varscan
	if( !defined( $h_cmd_log{"mpileup"}{$sample_id} ) ){
		&cmd_print("#&J	$sample_id\_mpileup	$sample_id\_aln");
		&cmd_print("$bin_samtools mpileup -AB -q 1 -d 50000 -Q 13 -f $ref_fa $sorted_bam -l $bed_region >$sorted_pileup");
		$h_cmd_log{"mpileup"}{$sample_id} = 1;
	}
	if( !defined( $h_cmd_log{"mpileup"}{$normal_sample_id} ) ){
		&cmd_print("#&J	$normal_sample_id\_mpileup	$normal_sample_id\_aln");
		&cmd_print("$bin_samtools mpileup -AB -q 1 -d 50000 -Q 13 -f $ref_fa $normal_sorted_bam -l $bed_region >$normal_sorted_pileup");
		$h_cmd_log{"mpileup"}{$normal_sample_id} = 1;
	}
	#varscan somatic
	my $last_job = "$sample_id\_mpileup,$normal_sample_id\_mpileup";
	&cmd_print("#&J	$sample_id\_VS_$normal_sample_id\_somaticVarscan	$last_job");
	&cmd_print("$bin_varscan somatic $normal_sorted_pileup $sorted_pileup $varscan_somatic_raw_vcf --output-vcf 1 --min-var-freq 0.0001 --somatic-p-value 0.2");

}


sub anno_all_variant{
	my $vcf_files_to_stat = "";
	my $anno_job_ids= "";
	#foreach my $sample_id( keys(%{$h_sample_info{"sample"}}) ){
		#foreach my $anno_type("GermlineGATK4","SomaticGATK4","SomaticVarDict"){
		#my $final_vcf = $h_sample_info{"call_variant"}{"GermlineGATK4"}{$sample_id}{"final_vcf"};
	foreach my $anno_type(sort {$a cmp $b} keys(%{$h_sample_info{"call_variant"}})){
		foreach my $sample_id(sort {$a cmp $b} keys(%{$h_sample_info{"call_variant"}{$anno_type}})){
			my $data_dir="$out_dir/smallVariant_anno/$anno_type/$sample_id";	`mkdir -p $data_dir`;
			my $to_anno_vcf = $h_sample_info{"call_variant"}{$anno_type}{$sample_id}{"final_vcf"};
			#print STDERR "$anno_type\t$sample_id\t$to_anno_vcf\n";
			my $left_align_vcf = "$data_dir/$sample_id.final.leftalign.vcf";
			my $yml_exomiser_file = "$data_dir/exomiser.$sample_id.yml";
			my $exomiser_anno_prefix ="$data_dir/exomiser.$sample_id";
			my $snpeff_anno_vcf="$data_dir/snpeff.$sample_id.vcf";
			my $snpeff_anno_stat="$data_dir/snpeff.$sample_id.stat.html";
			my $anno_gene_vcf="$data_dir/$sample_id.final.leftalign.annoGene.vcf";
			my $annoVar_out="$data_dir/$sample_id.final.leftalign.annoGene.annoVar";

			my $last_job = "$sample_id\_CallVariant$anno_type";
			&cmd_print("#&J	$sample_id\_AnnoVariant_$anno_type	$last_job");
			$anno_job_ids .= ",$sample_id\_AnnoVariant_$anno_type";
			&cmd_print("$bin_gatk LeftAlignAndTrimVariants --split-multi-allelics -R $ref_fa -V $to_anno_vcf -O $left_align_vcf");
			my $anno_exomiser_cmd = "";
			#if( defined( $h_sample_info{"call_variant"}{$anno_type}{$sample_id}{"hpo"} ) ){
			if( defined($h_sample_info{"sample"}{$sample_id}{"hpo"}) ){
				my $hpo_ids = join(",",keys(%{ $h_sample_info{"sample"}{$sample_id}{"hpo"} }));
				#my $hpo_ids = join(",",keys(%{$h_sample_info{"call_variant"}{$anno_type}{$sample_id}{"hpo"}}));
				&cmd_print("$bin_create_exomiser_exome liumm_vcf_file_path:$left_align_vcf liumm_output_prefix:$exomiser_anno_prefix liumm_hpo_ids:$hpo_ids >$yml_exomiser_file");
				&cmd_print("$bin_exomiser --analysis $yml_exomiser_file --spring.config.location=$config_exomiser");
				$anno_exomiser_cmd = "-exomiser $exomiser_anno_prefix.vcf";
			}
			&cmd_print("$bin_snpeff -s $snpeff_anno_stat $left_align_vcf >$snpeff_anno_vcf");
			&cmd_print("cat $snpeff_anno_vcf|$bin_fix_snpeff|$bin_select_snpeff|$bin_tag_nearbyOverlap_variant -ref $ref_fa|$bin_anno_gene $anno_exomiser_cmd -clinvar_all $database_clinvar_stat_all -clinvar_mul $database_clinvar_stat_mul -hpo_g2p $database_g2p_hpo -hpo_p2g $database_p2g_hpo >$anno_gene_vcf");
			#注释各类数据库信息by annovar
			&cmd_print("$bin_annovar_withPara -out $annoVar_out $anno_gene_vcf $annovar_humandb");
			$vcf_files_to_stat .= " $sample_id\_$anno_type:$anno_type:$annoVar_out\.hg19_multianno.vcf";
		}
	}
	my $data_dir = "$out_dir/final_result";	`mkdir -p $data_dir`;
	if( $anno_job_ids =~ s/^,// ){
		$vcf_files_to_stat =~ s/^,//;
		&cmd_print("#&J	AllStatSmallVairant	$anno_job_ids");
		&cmd_print("$bin_smallVariant $bed_region $data_dir/AllStat.smallVariant_stat.detail.xls $data_dir/AllStat.smallVariant_stat.stat.xls $vcf_files_to_stat");
	}
}


#stat sex hom VS het info
sub sex_stat{

	#my @samples = @_;
	if( defined( $h_cmd_log{"sex_stat"} ) ){
		print STDERR "ERR:only one sex_stat allowed\n";
		return;
	}
	else{
		$h_cmd_log{"sex_stat"} = 1;
	}

	my $last_job = "";
	my $input_vcfs = "";
	#foreach my $sample_id(@samples){
	foreach my $sample_id( keys(%{$h_sample_info{"sample"}}) ){
		if( defined($h_sample_info{"call_variant"}{"GermlineGATK4"}{$sample_id}{"final_vcf"}) ){
			$last_job .= ",$sample_id\_CallVariantGermlineGATK4";
			$input_vcfs .= " $sample_id:".$h_sample_info{"call_variant"}{"GermlineGATK4"}{$sample_id}{"final_vcf"};
		}
	}
	if($last_job =~ s/^,//){
		my $data_dir = "$out_dir/final_result";	`mkdir -p $data_dir`;
		&cmd_print("#&J	AllStatSex	$last_job");
		&cmd_print("$bin_statsex $input_vcfs >$data_dir/AllStat.sex_stat.xls");
	}
}

sub HLA_stat{

	if( defined( $h_cmd_log{"HLA_stat"} ) ){
		print STDERR "ERR:only one HLA_stat allowed\n";
		return;
	}
	else{
		$h_cmd_log{"HLA_stat"} = 1;
	}

	my $last_job = "";
	my $input_str = "";
	foreach my $sample_id( keys(%{$h_sample_info{"sample"}}) ){
		if( defined($h_sample_info{"sample"}{$sample_id}{"HLA"}) ){
			$last_job .= ",$sample_id\_HLA";
			$input_str .= " $sample_id:".$h_sample_info{"sample"}{$sample_id}{"HLA"};
		}
	}
	if($last_job =~ s/^,//){
		my $data_dir = "$out_dir/final_result";	`mkdir -p $data_dir`;
		&cmd_print("#&J	AllStatHLA	$last_job");
		&cmd_print("$bin_statHLA $data_dir/AllStat.HLA_stat.xls $input_str");
	}
}

sub msi_stat{

	if( defined( $h_cmd_log{"msi_stat"} ) ){
		print STDERR "ERR:only one msi_stat allowed\n";
		return;
	}
	else{
		$h_cmd_log{"msi_stat"} = 1;
	}

	my $last_job = "";
	my $input_msis = "";
	foreach my $sample_id( keys(%{$h_sample_info{"sample"}}) ){
		if( defined($h_sample_info{"sample"}{$sample_id}{"msi"}) ){
			$last_job .= ",$sample_id\_msisensor";
			$input_msis .= " $sample_id:".$h_sample_info{"sample"}{$sample_id}{"msi"};
		}
	}
	if($last_job =~ s/^,//){
		my $data_dir = "$out_dir/final_result";	`mkdir -p $data_dir`;
		&cmd_print("#&J	AllStatMSI	$last_job");
		&cmd_print("$bin_statmsi $data_dir/AllStat.msi_stat.xls $input_msis");
	}
}

sub factera_stat{

	if( defined( $h_cmd_log{"factera_stat"} ) ){
		print STDERR "ERR:only one factera_stat allowed\n";
		return;
	}
	else{
		$h_cmd_log{"factera_stat"} = 1;
	}

	my $last_job = "";
	my $input_str = "";
	foreach my $sample_id( keys(%{$h_sample_info{"sample"}}) ){
		if( defined($h_sample_info{"sample"}{$sample_id}{"factera_out"}) ){
			$last_job .= ",$sample_id\_factera";
			$input_str .= " $sample_id:".$h_sample_info{"sample"}{$sample_id}{"factera_out"};
		}
	}
	if($last_job =~ s/^,//){
		my $data_dir = "$out_dir/final_result";	`mkdir -p $data_dir`;
		&cmd_print("#&J	AllStatFactera	$last_job");
		&cmd_print("$bin_facterastat $input_str >$data_dir/AllStat.factera_stat.xls");
	}
}
