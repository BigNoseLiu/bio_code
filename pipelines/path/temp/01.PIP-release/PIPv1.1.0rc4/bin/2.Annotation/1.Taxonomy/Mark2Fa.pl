#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use FindBin qw($Bin);
use lib "$Bin/../../../lib/perl5";
use MyModule::GlobalVar qw($PIGZ $BLASTDB $BLASTN $SEQTAXID2TAX);

###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.0";
my $DATE = "2021-12-13";
###################################################################

my($inMark, $db, $outf, $fq, $outdir, $version);
GetOptions(
	"inMark:s" => \$inMark,
	"db:s" => \$db,
	"outf:s" => \$outf,
	"fq:s" => \$fq,
	"outdir:s" => \$outdir,
	"version|v" => \$version,
);
if($version){
	print basename $0." $VERSION\n";
	exit(0);
}
&help unless($inMark && $fq);
$db ||= "20210103";
$outf ||= "fasta";
$outdir ||= "./";
`mkdir -p $outdir` unless(-d $outdir);

open IN,"$SEQTAXID2TAX" or die $!;
my(%anno, %seqid2speid);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$seqid2speid{$t[0]} = $t[1];
	if($t[2] eq 'Viruses'){
		$anno{$t[0]} = $t[9];
	}else{
		$anno{$t[0]} = $t[8];
	}
}
close IN;

open IN,"$inMark" or die $!;
my %hash;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$t[1] =~ s/\/.*//g;
	$hash{$t[1]} = $t[2];
}
close IN;

open OU,">$outdir/validated.fna" or die $!;
foreach my $i(split /,/,$fq){
	my $name = basename $i;
	$name =~ s/\.gz//g; $name =~ s/\.fq|\.fastq//g;
	if($i =~ /\.fq\.gz$/ || $i =~ /\.fastq\.gz$/){
        open FQ,"$PIGZ -dc -p 16 $i |" or die $!;
    }elsif($i =~ /\.fq$/ || $i =~ /\.fastq$/){
        open FQ,"$i" or die $!;
    }else{
        die "Please confirm the format of input reads file!\n";
    }
	open OF,"| $PIGZ -c -p 16 >$outdir/$name.Anno.fna.gz" or die $!;
	while(<FQ>){
		chomp;
		my $id = $_;
		my $taxid = (split /\|/,$id)[-1];
		chomp(my $seq = <FQ>);
		<FQ>;
		chomp(my $qual = <FQ>);
		$id =~ s/^@//g;
		my $id2 = $id;
		$id2 =~ s/\/.*//g;
		$id2 =~ s/ .*//g;
		$id =~ s/ /|/g;
		if(exists $hash{$id2}){
			if($outf eq 'fasta'){
				print OU ">$id\n$seq\n";
			}else{
				print OU "\@${id}\n$seq\n+\n$qual\n";
			}
		}
		if(exists $anno{$taxid}){
			print OF ">${id} $anno{$taxid}\n$seq\n";
		}elsif($taxid == 77643){
			print OF ">${id} Mycobacterium tuberculosis complex\n$seq\n";
		}
	}
	close FQ;
	close OF;
}
close OU;

`$BLASTN -query $outdir/validated.fna -db $BLASTDB/$db/nt -num_threads 20 -out $outdir/validated.m8 -outfmt \"6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp\" -evalue 0.05 -word_size 28 -reward 1 -penalty -2 -perc_identity 95 -qcov_hsp_perc 90`; #-num_alignments 100000

open IN,"$outdir/validated.m8" or die $!;
my %score;
while(<IN>){
    chomp;
    my @t = split /\t/,$_;
    if(!exists $score{$t[0]} || $score{$t[0]}{"identity+coverage"} < $t[2]+$t[14]){
        $score{$t[0]}{"identity+coverage"} = $t[2] + $t[14];
    }
}
close IN;

open IN,"$outdir/validated.m8" or die $!;
my %mark;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	my $queryTaxID = (split /\|/,$t[0])[-1];
	my $targetTaxID;
	if($t[1] =~ /:/){
		$targetTaxID = (split /:/,$t[1])[0];
	}else{
		$targetTaxID = (split /\|/,$t[1])[0];
	}
	next unless($targetTaxID);
	if($t[2] + $t[14] == $score{$t[0]}{"identity+coverage"}){
		if(exists $seqid2speid{$queryTaxID} && exists $seqid2speid{$targetTaxID} && $seqid2speid{$queryTaxID} == $seqid2speid{$targetTaxID}){
			$mark{$t[0]} = 1;
		}
	}
}
close IN;

open VO,">$outdir/validated.out" or die $!;
foreach my $i(keys %score){
	my $queryTaxID = (split /\|/,$i)[-1];
	if(exists $mark{$i}){
		print VO "$i\t$anno{$queryTaxID}\ttrue\n" if(exists $anno{$queryTaxID});
	}else{
		print VO "$i\t$anno{$queryTaxID}\tfalse\n" if(exists $anno{$queryTaxID});
	}
}

# Fix some missing verification results
open IN,"$outdir/validated.fna" or die $!;
while(<IN>){
	chomp;
	next unless(/^>/);
	my $id = $_;
	$id =~ s/^>//;
	my $queryTaxID = (split /\|/,$id)[-1];
	print VO "$id\t$anno{$queryTaxID}\tfalse\n" if(!exists $score{$id} && exists $anno{$queryTaxID});
}
close IN;
close VO;

sub help{
print "
		Usage: perl $0

		--inMark       Enter a list of sequences to be validated.
		--db           Blast database version. [20210103]
		--outf         Output file's format, [fasta]|fastq
		--fq           Fastq file, separated with comma, example: --fq fq1,fq2
		--outdir       Output path [./]
		--version|v    Print version information 
\n";
exit(0);
}
