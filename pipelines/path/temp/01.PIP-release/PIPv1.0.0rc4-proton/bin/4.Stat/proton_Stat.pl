#!/usr/bin/env perl

#===============================================================================
#   FileName: Stat.pl
#   Author  : qiuwancen
#   E-Mail  : 972538446@qq.com
#   Version : 1.0
#   Date    : 2022-03-22
#   Description: Copy the result files to Result directory.
#===============================================================================

use warnings;
use strict;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use lib "$Bin/../../lib/perl5";
use MyModule::MyLog qw(printInfo2ERR);
use MyModule::GlobalVar qw($TXT2EXCEL $PERL $PYTHON3 $KRAKENDB);


# hy :增加$seqtype
my($list, $analysis, $workdir, $result, $shellDir, $help, $seqtype);
GetOptions(
    "analysis=s" => \$analysis,
    "workdir=s"  => \$workdir,
    "result=s"   => \$result,
    "shellDir=s" => \$shellDir,
    "type=s" => \$seqtype,
    "help|h"     => \$help,
);
pod2usage(-exitval => 1, -verbose => 1) if($help || @ARGV < 1);
($list, ) = @ARGV;
$list = abs_path $list;


#===============================================================================
#   Set the default values of options
#===============================================================================
$analysis ||= "Merge-QC:Fastp,RemoveHost,rRNAFilter-Annotation:Taxonomy,ResistanceGene,VFDB";
$seqtype ||= "SE";
$workdir  ||= "./";
$shellDir ||= "./";
$workdir  = abs_path $workdir;
$result   =   "$workdir/plugin_out/Result";
`mkdir -p $result` unless(-d $result);
$result   = abs_path $result;
$shellDir = abs_path $shellDir;

#===============================================================================
#   Create the analytic directory
#===============================================================================
#`mkdir -p $shellDir` unless(-d $shellDir);
`mkdir -p $workdir/basecaller_results` unless(-d "$workdir/basecaller_results");

#===============================================================================
#   Scripts
#===============================================================================
my $QC            = "$Bin/proton_QC.stat.pl"; #hy
my $Normalize     = "$Bin/Normalize.py";
my $CommonSpecies = "$Bin/CommonSpecies.pl";
my $mergeResults  = "$Bin/proton_mergeResults.py"; #hy
my $fillCoverage  = "$Bin/fillCoverage.py";
my $plot          = "$Bin/proton_Plot.py"; #hy
my $split_qc      = "$Bin/split_QC_Stat.py";

#===============================================================================
#   Global variables
#===============================================================================
## Directory
my %path = (
    "QC"                => "$workdir/01.QC",
    "Taxonomy"          => "$workdir/02.Annotation/01.Taxonomy",
    "ResistanceGene"    => "$workdir/02.Annotation/02.ResistanceGene",
    "VFDB"              => "$workdir/02.Annotation/03.VFDB",
);

#===============================================================================
#   Get the sample list
#===============================================================================

#hy -start 20220511
my %config;
parse_config("$Bin/../../conf/configure.ini",\%config);
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		if (/(\S+)=(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;
		}
	}
	close IN;
}
#hy -end

my %samp_list = ();
open IN,"$list" or die $!;
my $flag = 0;
while (<IN>) {
    chomp;
	if($_=~/^\s*$/){next;} #hy
	my @c=split /\t/,$_; #hy
    my($name, $type, $fq) = (split /\t/,$_)[0, 1, 2];
	if($flag == 0){
		my $dir = dirname $fq;
		#`cp $dir/*.summaryReport.html $workdir/basecaller_results`; #hy
		
		#hy -start
		if($seqtype eq "SE"){
			`$Bin/FastqCount-master/FastqCount_v0.5 $c[2] >$workdir/00.Merge/$name/$type/sta.log`; #hy 20220511
		}elsif($seqtype eq "PE"){
			`$config{pigz} -dc $c[2] $c[3] | $Bin/FastqCount-master/FastqCount_v0.5 - >$workdir/00.Merge/$name/$type/sta.log`; #hy 20220511
		}
		#hy -end
		
	}
    my $sample = "$name\_$type";
    $samp_list{$name}{$type} = $sample;
}
close IN;


#===============================================================================
#   Generate the scripts
#===============================================================================
`$PERL $QC $path{QC}/sample.list $result/QC.Stat`;
`$PERL $TXT2EXCEL $result/QC.Stat $result/QC.Stat.xlsx`;
`$PYTHON3 $Normalize --qc $result/QC.Stat --anno $path{Taxonomy}/sample.list`;
`$PERL $CommonSpecies --list $path{Taxonomy}/sample.list --outdir $result`;
#`$PYTHON3 $plot --qc $result/QC.Stat --profile $result/Pathogeny_CommonSpecies.Normalize.txt --outdir $result`;
chdir("$result");
`zip -r $result/zhikong.zip Total_CommonSpecies* Pathogeny_CommonSpecies*`;

for my $name (keys %samp_list) {
    my $sampDir = "$result/$name";
    `mkdir -p $sampDir` unless(-d $sampDir);
	`rm $workdir/by2.sample` if(-e "$workdir/by2.sample");

    # Copy the files of identification
    if ($analysis =~ /Taxonomy/) {
    	my @tmp = split(/=/, $analysis);
    	my $db  = (split(/,/, $tmp[1]))[0];
        my $origTaxDir = $path{Taxonomy};
        for my $type (keys %{$samp_list{$name}}) {
			`$PYTHON3 $fillCoverage --db $KRAKENDB/$db --indir $workdir --name $name --library $type --outdir $sampDir`;
			`cp -fr $origTaxDir/$name/$type/*.xls $sampDir`;
			`cp -fr $origTaxDir/$name/$type/Pathogeny_Detail.xlsx $sampDir`;
			`cp -fr $origTaxDir/$name/$type/Pathogeny_Detail.txt $sampDir`;
			`cp -fr $origTaxDir/$name/$type/Total_Detail.xlsx $sampDir`;
			`cp -fr $origTaxDir/$name/$type/Total_Detail.txt $sampDir`;
			`cp -fr $origTaxDir/$name/$type/Taxonomy_Summary.txt $sampDir`;
			`cp -fr $origTaxDir/$name/$type/bacteria.result.fq.gz $sampDir/`; ###temp
			`cp -fr $origTaxDir/$name/$type/fungi.result.fq.gz $sampDir/`; ###temp
			`cp -fr $origTaxDir/$name/$type/protozoa.result.fq.gz $sampDir/`; ###temp
			`cp -fr $origTaxDir/$name/$type/viral.result.fq.gz $sampDir/`; ###temp
			`cp -fr $origTaxDir/$name/$type/LY.result.fq.gz $sampDir/`; ###temp
			`cp -fr $sampDir/bacteria.xls $sampDir/bacteria.pick.xls`;
			`cp -fr $sampDir/fungi.xls $sampDir/fungi.pick.xls`;
			`cp -fr $sampDir/LY.xls $sampDir/LY.pick.xls`;
			`cp -fr $sampDir/protozoa.xls $sampDir/protozoa.pick.xls`;
			`cp -fr $sampDir/viral.xls $sampDir/viral.pick.xls`;
			`$PYTHON3 $mergeResults --indir $workdir --name $name --library $type --outdir $sampDir`;
        }
    }
    
    # Copy the files of VF identification
    if ($analysis =~ /VFDB/) {
        my $origVfDir = $path{"VFDB"};
        for my $type (keys %{$samp_list{$name}}) {
            `cp -fr $origVfDir/$name/$type/VF.xls $sampDir`;
        }
    }

	# Copy the files of ResistanceGene identification
	if ($analysis =~ /ResistanceGene/) {
		my $origResGeneDir = $path{"ResistanceGene"};
		for my $type (keys %{$samp_list{$name}}){
			`cp -fr $origResGeneDir/$name/$type/ARG.cov.pick.xls $sampDir`;
		}
	}
    
    # Copy the files of RNA vs DNA
#    if($analysis =~ /RNADNA/){
#        my $origVsDir = $path{"RNADNA"};
#        my $destVsDir = $destPath{"RNADNA"};
#        system("mkdir -p $destVsDir") unless (-e $destVsDir);
#        print SH "cp -fr $origVsDir/Result/$name/* $destVsDir\n";
#    }

	# Copy the files of SpeciesExtract
#	if($analysis =~ /SpeciesExtract/){
#		my $origExtractDir = $path{"SpeciesExtract"};
#		my $destExtractDir = $destPath{"SpeciesExtract"};
#		for my $type (keys %{$samp_list{$name}}){
#			my $destTypeDir = "$destExtractDir/$type";
#			if(-d "$origExtractDir/$name/$type"){
#				system("mkdir -p $destTypeDir") unless (-e $destTypeDir);
#				print SH "cp -fr $origExtractDir/$name/$type/* $destTypeDir\n";
#			}
#		}
#	}
}
#close SH;

###############################################################################################################
foreach my $smp(keys %samp_list){
	`cp $result/$smp/Total_Detail.xlsx $result/$smp/$smp.merge.xls`
}

`$PYTHON3 $split_qc -d $workdir`;

open RUN,">$workdir/by2.run" or die $!;
print RUN "Result\n";
for my $name (keys %samp_list) {
	print RUN "$name\n";
}
close RUN;

`sleep 30`;

`touch $workdir/by2.sample`;
`touch $workdir/by2.done`;

#===============================================================================
#   Execute the scripts
#===============================================================================
#printInfo2ERR("Start program $0");
#`sh $shellDir/upload.sh`;
#printInfo2ERR("Finish program $0");

__END__

=pod

=head1 NAME

upload.pl - The program for uploading result

=head1 VERSION

v1.0

=head1 SYNOPSIS

perl upload.pl <sample.list> [options]

=head1 ARGUMENTS

=over 8

=item B<sample.list> <file>

List of raw Reads.

Format: "sample DNA|RNA read1.fq [read2.fq]"

There are 4 columns separated with Tab in the file, the 1st is sample name, the
2nd is sample type (DNA or RNA), the 3rd is the fastq file of reads (read1 for PE),
the 4th is the fastq file of read2 (only for PE).

=back

=head1 OPTIONS

=over 8

=item B<--analysis> <Str>

The analytic steps you done. The different steps separated by '-', the
arguments for every step follow the ':' and separated by ','. The selectable
steps as bollow:
default: [Merge-QC:Fastp,RemoveHost,rRNAFilter-Annotation:Taxonomy,ResistanceGene,VFDB]

=item B<--workdir> <Path>

The path including 00.Rawdata,01.QC,.... [./]

=item B<--result> <Path>

Result path. [outdir/Result]

=item B<--shellDir> <path>

The directory to store the scripts generated by this program. [./]

=item B<--help|h>

Print this information.

=back

=cut
