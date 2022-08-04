#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename;

###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.0";
my $DATE = "2021-03-28";
###################################################################

my($ko, $fqs, $outdir, $version);
GetOptions(
	'krakenout:s' => \$ko,
	'fq:s' => \$fqs,
	'outdir:s' => \$outdir,
	'version|v' => \$version,
);

if($version){
    print basename $0." $VERSION\n";
    exit(0);
}

&help unless($ko && $fqs);
$outdir ||= "./";
`mkdir -p $outdir` unless(-d $outdir);
$outdir = abs_path $outdir;

my @FQs = split /,/, $fqs;
if(scalar @FQs == 1){
	if($FQs[0] =~ /\.fq\.gz$/ || $FQs[0] =~ /\.fastq\.gz$/){
		open FQ, "gzip -dc $FQs[0] |" or die $!;
	}elsif($FQs[0] =~ /\.fq$/ || $FQs[0] =~ /\.fastq$/){
		open FQ, "$FQs[0]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	my $name = basename $FQs[0];
	$name =~ s/\.gz$//g;
	open OU,"| gzip > $outdir/$name.gz" or die $!;
	open KO,"$ko" or die $!;
	my $id2 = 0;
	my($FqID, $Seq, $Qual);
	while(<KO>){
		chomp;
		next unless(/^U/);
		my $Rid = (split /\t/,$_)[1];
		$Rid =~ s/\/.*//g;
		while($id2 ne $Rid){
			chomp($FqID = <FQ>);
			$id2 = $FqID;
			$id2 =~ s/ .*//g; $id2 =~ s/\/.*//g; $id2 =~ s/^\@//g;
			chomp($Seq = <FQ>);
			<FQ>;
			chomp($Qual = <FQ>);
		}
		print OU "$FqID\n$Seq\n+\n$Qual\n";
	}
	close KO; close OU;
}else{
	if($FQs[0] =~ /\.fq\.gz$/ || $FQs[0] =~ /\.fastq\.gz$/){
		open FQ1,"gzip -dc $FQs[0] |" or die $!;
	}elsif($FQs[0] =~ /\.fq$/ || $FQs[0] =~ /\.fastq$/){
		open FQ1,"$FQs[0]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	if($FQs[1] =~ /\.fq\.gz$/ || $FQs[1] =~ /\.fastq\.gz$/){
		open FQ2,"gzip -dc $FQs[1] |" or die $!;
	}elsif($FQs[1] =~ /\.fq$/ || $FQs[1] =~ /\.fastq$/){
		open FQ2,"$FQs[1]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	my $name1 = basename $FQs[0];   $name1 =~ s/\.gz//g;
	my $name2 = basename $FQs[1];   $name2 =~ s/\.gz//g;
	open OU1,"| gzip > $outdir/$name1.gz" or die $!;
	open OU2,"| gzip > $outdir/$name2.gz" or die $!;
	open KO,"$ko" or die $!;
	my($FqID1, $FqID2, $FqID2_2, $Seq1, $Seq2, $Qual1, $Qual2);
	my $FqID1_2 = 0;
	while(<KO>){
		chomp;
		next unless(/^U/);
		my $Rid = (split /\t/,$_)[1];
		$Rid =~ s/\/.*//g;
		while($FqID1_2 ne $Rid){
			chomp($FqID1 = <FQ1>);    chomp($FqID2 = <FQ2>);
			chomp($Seq1 = <FQ1>);     chomp($Seq2 = <FQ2>);
			<FQ1>;                    <FQ2>;
			chomp($Qual1 = <FQ1>);    chomp($Qual2 = <FQ2>);
			$FqID1_2 = $FqID1;
			$FqID1_2 =~ s/ .*//g; $FqID1_2 =~ s/\/.*//g; $FqID1_2 =~ s/^\@//g;
		}
		print OU1 "$FqID1\n$Seq1\n+\n$Qual1\n";
		print OU2 "$FqID2\n$Seq2\n+\n$Qual2\n";
	}
	close KO; close OU1; close OU2;
}

sub help{
print "
		Usage: perl $0

		--krakenout    intermediate file of kraken analysis, the file's suffix is kout
		--fq           fastq file, separated with comma, example: --fq fq1,fq2
		--outdir       output path [.]
		--version|v    print version information.

";
exit(0)
}
