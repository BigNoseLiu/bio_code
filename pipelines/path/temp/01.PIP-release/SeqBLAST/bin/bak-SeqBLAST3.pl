#!/usr/bin/env perl

#===============================================================================
#   Author  : qiuwancen
#   E-Mail  : 972538446@qq.com
#   Version : 1.0
#   Date    : 2022-05-12
#===============================================================================

use warnings;
use strict;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use lib "$Bin/../lib";
use MyModule::GlobalVar qw($PIGZ $BLASTN $BLASTDB);

my($species, $fq, $dbVersion, $module, $num, $outdir);
GetOptions(
	"species:s" => \$species,
	"fq:s" => \$fq,
	"dbVersion:s" => \$dbVersion,
	"module:s" => \$module,
	"num:s" => \$num,
	"outdir:s" => \$outdir,
);
&help unless($species && $fq);
$dbVersion ||= "20210103";
$module ||= "BLAST";
$num ||= 20;
$outdir ||= "./";
`mkdir -p $outdir` unless(-d $outdir);
################################################################################
if($module eq 'BLAST'){
	open IN,"$PIGZ -dc $fq |" or die $!;
	my $spe = $species;
#	$spe =~ s/[\/ ]/_/g; $spe =~ s/[\[()\]]//g;
	$spe =~ s/[\/?#\$%^&*@!~|\\()\[\] '"]/_/g;
	open OU,">$outdir/$spe.fna" or die $!;
	while(<IN>){
		chomp;
		last if($num == 0);
		my $id = $_;
		chomp(my $seq = <IN>);
		if($id =~ /$species$/ && $num > 0){
			print OU "$id\n$seq\n";
			$num--;
		}
	}
	close IN;
	close OU;

	`$BLASTN -query $outdir/$spe.fna -db $BLASTDB/$dbVersion/nt -out $outdir/$spe.m8 -outfmt "6 qseqid sacc pident length evalue bitscore qcovhsp" -num_threads 12 -evalue 0.05 -word_size 28 -reward 1 -penalty -2`;
	`$Bin/m82annoTaxid.pl --inM8 $outdir/$spe.m8 --out $outdir/$spe.m8.anno`;
	`echo "All done." > $outdir/blast.done`;
}elsif($module eq 'EXTRACT'){
	open IN,"$PIGZ -dc $fq |" or die $!;
	my $spe = $species;
#	$spe =~ s/[\/ ]/_/g; $spe =~ s/[\[()\]]//g;
	$spe =~ s/[\/?#\$%^&*@!~|\\()\[\] '"]/_/g;
	open OU,">$outdir/$spe.fna" or die $!;
	while(<IN>){
		chomp;
		my $id = $_;
		chomp(my $seq = <IN>);
		if($id =~ /$species$/){
			print OU "$id\n$seq\n";
		}
	}
	close IN;
	close OU;
	`echo "All done." > $outdir/extract.done`;
}

sub help{
print "
		Usage: perl $0

		--species      species name
		--fq           input fastq file
		--dbVersion    NT version
		--module       [BLAST]|EXTRACT
		--num          extract the first [INT, 20] sequences
		--outdir       [./]

";
exit(0);
}
