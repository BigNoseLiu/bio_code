#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use FindBin qw($Bin);
use lib "$Bin/../../../lib/perl5";
use MyModule::GlobalVar qw($BIN_PATH $SRC);

###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.0";
my $DATE = "2021-09-16";
###################################################################

my($report, $db, $level, $fq, $krakenOut, $outdir, $version);
GetOptions(
	"report:s" => \$report,
	"db:s" => \$db,
	"level:s" => \$level,
	"fq:s" => \$fq,
	"krakenout:s" => \$krakenOut,
	"outdir:s" => \$outdir,
	"version|v" => \$version,
);
if($version){
    print basename $0." $VERSION\n";
    exit(0);
}
&help unless($report && $db && $fq && $krakenOut);
$level  ||= "Family";
$outdir ||= "./";
$outdir = abs_path $outdir;
`mkdir -p $outdir` unless(-d $outdir);

##########################################################
my $getSpeciesReads = "$Bin/getSpeciesReads.pl";
##########################################################

open IN,"$db/SeqTaxid2${level}Taxid-K.patho.list" or die $!;
my %taxid;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$taxid{$t[1]}{kingdom} = $t[0];
	$taxid{$t[1]}{$level} = $t[2];
}
close IN;

open RP,"$report" or die $!;
my %save;
while(<RP>){
	chomp;
	my @t = split /\t/,$_;
	if(exists $taxid{$t[4]}){
		push @{$save{$taxid{$t[4]}{kingdom}}{$taxid{$t[4]}{$level}}{id}}, $t[4];
	}
}
close RP;

open OU,">$outdir/$level.taxid.list" or die $!;
foreach my $k(sort keys %save){
	foreach my $g(keys %{$save{$k}}){
		print OU "$k\t$g\t".(join ",", @{$save{$k}{$g}{id}})."\n";
	}
}
close OU;

`$getSpeciesReads --fq $fq --taxid $outdir/$level.taxid.list --krakenout $krakenOut --db $db --outf fasta --outdir $outdir/Fasta`;
`$SRC/watchDog.pl --mem 2G --num_paral 2 $outdir/Fasta/Align.sh`;
#`find $outdir/Fastq -name "unmapped" |xargs cat >> $outdir/Fastq/delete.list`;
#`sed -i 's/^-e //g' $outdir/Fastq/delete.list`;

#open IN,"$outdir/Fastq/delete.list";
#my %delete;
#while(<IN>){
#	chomp;
#	my @t = split /\t/,$_;
#	if(exists $save{$t[0]}{$t[1]}){
#		foreach my $i(@{$save{$t[0]}{$t[1]}{id}}){
#			$delete{$i} = 1;
#		}
#	}
#}
#close IN;

#open RP,"$report" or die $!;
#open OU,">$report.2" or die $!;
#while(<RP>){
#	chomp;
#	my @t = split /\t/,$_;
#	if(!exists $delete{$t[4]}){
#		print OU "$_\n";
#	}
#}
#close RP;
#close OU;
#`mv $report.2 $report`;

sub help{
print "
	Usage: perl $0

	--report             kraken's report file
	--db                 path of the database to annotate
	--fq                 clean data's fastq file, separated with comma, example: --fq fq1,fq2
	--krakenout          intermediate file of kraken analysis, the file's suffix is out
	--outdir             output path. [./]
	--version|v          print version information.

";
exit(0);
}
