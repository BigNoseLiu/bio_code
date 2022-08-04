#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Basename;
use FileHandle;
use lib "$Bin/../../../lib/perl5";
use MyModule::GlobalVar qw($PIGZ);

###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.0";
my $DATE = "2021-02-02";
###################################################################
my($fqs, $outf, $taxid, $krakenOut, $outdir, $version);
GetOptions(
	"fq:s" => \$fqs,
	"outf:s" => \$outf,
	"taxid:s" => \$taxid,
	"krakenout:s" => \$krakenOut,
	"outdir:s" => \$outdir,
	"version|v" => \$version,
);
if($version){
    print basename $0." $VERSION\n";
    exit(0);
}
&help unless($fqs && $taxid && $krakenOut);
$outf ||= "fastq";
$outdir ||= "./Result";
`mkdir -p $outdir` unless(-d $outdir);
$outdir = abs_path $outdir;

open TAX,"$taxid" or die $!;
my(%s2k, %id2s, %fh, @S);
my @FQs = split /,/, $fqs;
while(<TAX>){
	next unless(/^[BEAV]/);
	chomp;
	my @t = split /\t/,$_;
	$t[1] =~ s/ /_/g;	$t[1] =~ s/[\(\)']//g;
	$s2k{$t[1]} = $t[0];
	foreach my $i(split /,/, $t[2]){
		$id2s{$i} = $t[1];
	}
	push @S, $t[1];
}
close TAX;

if(scalar @FQs == 1){
	if($FQs[0] =~ /\.fq\.gz$/ || $FQs[0] =~ /\.fastq\.gz$/){
		open FQ,"$PIGZ -dc -p 16 $FQs[0] |" or die $!;
	}elsif($FQs[0] =~ /\.fq$/ || $FQs[0] =~ /\.fastq$/){
		open FQ,"$FQs[0]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	my $name = basename $FQs[0];	$name =~ s/\.gz//g;
	if($outf eq 'fastq'){
		for(@S){
			`mkdir -p "$outdir/$s2k{$_}/$_"` unless(-d "$outdir/$s2k{$_}/$_");
			open $fh{$_},"| $PIGZ -c -p 16 >$outdir/$s2k{$_}/$_/$name.gz" or die $!;
		}
	}elsif($outf eq 'fasta'){
		$name =~ s/\.fq|\.fastq//g;
		for(@S){
			`mkdir -p "$outdir/$s2k{$_}/$_"` unless(-d "$outdir/$s2k{$_}/$_");
			open $fh{$_},"| $PIGZ -c -p 16 >$outdir/$s2k{$_}/$_/$name.fa.gz" or die $!;
		}
	}else{
		die "Please confirm the format of output file!\n";
	}

	open KO,"$krakenOut" or die $!;
	my($fqID, $fqSeq, $fqQual);
	my $fqID2 = 0;
	while(<KO>){
		chomp;
		my($Rid, $Tid) = (split /\t/,$_)[1,2];
		$Rid =~ s/\/.*//g;
		next unless(exists $id2s{$Tid});
		while($fqID2 ne $Rid){
			chomp($fqID = <FQ>);
			chomp($fqSeq = <FQ>);
			<FQ>;
			chomp($fqQual = <FQ>);
			$fqID =~ s/^\@//g;
			$fqID2 = $fqID;
			$fqID2 =~ s/\/.*$//g;	$fqID2 =~ s/ .*//g;
		}
		if($outf eq 'fastq'){
			$fh{$id2s{$Tid}}->print("\@$fqID\n$fqSeq\n+\n$fqQual\n");
		}else{
			$fh{$id2s{$Tid}}->print(">$fqID\n$fqSeq\n");
		}
	}
	close KO;

	for(keys %s2k){
		close $fh{$_};
	}
}else{
	if($FQs[0] =~ /\.fq\.gz$/ || $FQs[0] =~ /\.fastq\.gz$/){
		open FQ1,"$PIGZ -dc -p 16 $FQs[0] |" or die $!;
	}elsif($FQs[0] =~ /\.fq$/ || $FQs[0] =~ /\.fastq$/){
		open FQ1,"$FQs[0]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	if($FQs[1] =~ /\.fq\.gz$/ || $FQs[1] =~ /\.fastq\.gz$/){
		open FQ2,"$PIGZ -dc -p 16 $FQs[1] |" or die $!;
	}elsif($FQs[1] =~ /\.fq$/ || $FQs[1] =~ /\.fastq$/){
		open FQ2,"$FQs[1]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	my $name1 = basename $FQs[0];	$name1 =~ s/\.gz//g;
	my $name2 = basename $FQs[1];	$name2 =~ s/\.gz//g;
	if($outf eq 'fastq'){
		for(@S){
			open $fh{$_}{1},"| $PIGZ -c -p 16 >$outdir/$s2k{$_}/$_/$name1.gz" or die $!;
			open $fh{$_}{2},"| $PIGZ -c -p 16 >$outdir/$s2k{$_}/$_/$name2.gz" or die $!;
		}
	}elsif($outf eq 'fasta'){
		$name1 =~ s/\.fq|\.fastq//g;	$name2 =~ s/\.fq|\.fastq//g;
		for(@S){
			open $fh{$_}{1},"| $PIGZ -c -p 16 >$outdir/$s2k{$_}/$_/$name1.fa.gz" or die $!;
			open $fh{$_}{2},"| $PIGZ -c -p 16 >$outdir/$s2k{$_}/$_/$name2.fa.gz" or die $!;
		}
	}else{
		die "Please confirm the format of output file!\n";
	}
	
	open KO,"$krakenOut" or die $!;
	my($fqID1, $fqID2, $fqID2_2, $fqSeq1, $fqSeq2, $fqQual1, $fqQual2);
	my $fqID1_2 = 0;
	while(<KO>){
		chomp;
		my($Rid, $Tid) = (split /\t/,$_)[1,2];
		next unless(exists $id2s{$Tid});
		while($fqID1_2 ne $Rid){
			chomp($fqID1 = <FQ1>);    chomp($fqID2 = <FQ2>);
			chomp($fqSeq1 = <FQ1>);   chomp($fqSeq2 = <FQ2>);
			<FQ1>;                    <FQ2>;
			chomp($fqQual1 = <FQ1>);  chomp($fqQual2 = <FQ2>);
			$fqID1 =~ s/^\@//g;
			$fqID1_2 = $fqID1;	$fqID1_2 =~ s/\/.*$//g;   $fqID1_2 =~ s/ .*//g;
		}
		if($outf eq 'fastq'){
			$fh{$id2s{$Tid}}{1}->print("\@$fqID1\n$fqSeq1\n+\n$fqQual1\n");
			$fh{$id2s{$Tid}}{2}->print("\@$fqID2\n$fqSeq2\n+\n$fqQual2\n");
		}else{
			$fh{$id2s{$Tid}}{1}->print(">$fqID1\n$fqSeq1\n");
			$fh{$id2s{$Tid}}{2}->print(">$fqID2\n$fqSeq2\n");
		}
	}
	close KO;

	for(keys %s2k){
		close $fh{$_}{1};
		close $fh{$_}{2};
	}
}

sub help{
print "
		Usage: perl $0

		Mandatory Parameters:
			--fq             clean data's fastq file, separated with comma, example: --fq fq1,fq2
			--taxid          file that including a list of species taxid separated with comma,
                             example: id1,id2,...
			--krakenout      intermediate file of kraken analysis, the file's suffix is out

		Optional Parameters:
			--outf           output file's format, [fastq]|fasta
			--outdir         output path [./Result]
			--version|v      print version information.

		Example:
			perl $0 --fq fq1,fq2 --taxid id1,id2,id3 --krakenout sample.out --outdir Result
			perl $0 --fq fq1 --outf fasta --taxid id --krakenout sample.out

";
exit(0);
}
