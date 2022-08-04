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
my $DATE = "2022-03-08";
###################################################################

my($kout, $top, $db, $out, $kmerCut);
GetOptions(
	"kout:s" => \$kout,
	"top:s" => \$top,
	"db:s" => \$db,
	"out:s" => \$out,
	"kmerCut:s" => \$kmerCut,
);
&help unless($kout && $db);
$top ||= 5;
$out ||= "./out.txt";
$kmerCut ||= 0.6875;

open IN,"$db/all.Taxonomy.txt" or die $!;
my %patho;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	# Correlate the taxID of the species with kingdom's information  
	$patho{$t[1]} = $t[2] if($t[12] eq '*');
}
close IN;

open IN,"$db/seqTaxid2speciesTaxid" or die $!;
my %seqid2speid;
while(<IN>){
	next if(/NULL/);
	chomp;
	my @t = split /\t/,$_;
	$seqid2speid{$t[0]} = $t[1];
	if(exists $patho{$t[1]}){
		$patho{$t[0]} = $patho{$t[1]};
	}
}
close IN;

###################################
# 获取seqtaxid对应的S1水平taxid、S水平taxid
# 表头：Seqtid\tS1_Taxid\tS_Taxid
open IN,"$db/seqTaxid2S1LevelTaxid" or die $!;
my (%seqtid2S1taxid, %seqtid2Stid);
while(<IN>){
	next if(/NULL/);
	chomp;
	my @t = split /\t/,$_;
	$seqtid2S1taxid{$t[0]} = $t[1];
	$seqtid2Stid{$t[0]} = $t[2];
	if(exists $patho{$t[2]}){
		$patho{$t[0]} = $patho{$t[2]};
	}
}
close IN;

##########打乱文件行的顺序-消除模拟数据的局限（模拟数据 同一条染色体的reads的命名相似，容易在文件中排在一起） ---huanghy
my $backupfile="$kout"."_backup";
if(! -s $backupfile){
	`mv $kout $backupfile`;
}
`shuf $backupfile -o $kout`;  #打乱文件行顺序
##################################################


open IN,"$kout" or die $!;
open OU,">$out" or die $!;
my %sort;
while(<IN>){
	next if(/^U/);
	chomp;
	my @t = split /\t/,$_;
	#next if($patho{$t[2]} eq 'Viruses');
	next if($t[2] == 9606);
	if(exists $seqid2speid{$t[2]}){
		my($IDnum, $TotalKmer);
		my %uniqspe;
		foreach my $i(split / /,$t[4]){
			my($taxid, $num) = split /:/,$i;
			# calculate the total kmers count of one read 
			$TotalKmer += $num;
			if($patho{$t[2]} eq 'Viruses' and exists $seqtid2S1taxid{$t[2]}){
				$IDnum += $num if($seqtid2S1taxid{$taxid} == $seqtid2S1taxid{$t[2]} || $taxid == $seqid2speid{$t[2]});
				if(exists $seqtid2S1taxid{$taxid}){
					$uniqspe{$seqtid2S1taxid{$taxid}} = 1;
				}
			}else{
				$IDnum += $num if($seqid2speid{$taxid} == $seqid2speid{$t[2]});
				if(exists $seqid2speid{$taxid}){
					$uniqspe{$seqid2speid{$taxid}} = 1;
				}
			}
		}	    
		# Outputs a unique comparison of species
		if(scalar(keys %uniqspe) == 1){
			my $Rate = $IDnum / $TotalKmer;
			print OU "$_\t$IDnum\t$Rate\n";
			# Count the number of corresponding Kmer confidence scores
			if(exists $seqtid2S1taxid{$t[2]}){
				$sort{$seqtid2S1taxid{$t[2]}}{$Rate} += 1;
			}else{
				$sort{$seqid2speid{$t[2]}}{$Rate} += 1;
			}
		}elsif(exists $seqtid2S1taxid{$t[2]}){
			my $Rate = $IDnum / $TotalKmer;
			print OU "$_\t$IDnum\t$Rate\n";
			if(exists $seqtid2S1taxid{$t[2]}){
				$sort{$seqtid2S1taxid{$t[2]}}{$Rate} += 1;
			}else{
				$sort{$seqid2speid{$t[2]}}{$Rate} += 1;
			}
		}
	}
}
close IN;
close OU;

# Calculate the minimum confidence score of the sequence with the top N confidence score
my $outdir = dirname(abs_path $out);
open NC,">$outdir/Species-HighConfidenceKmerScore.txt" or die $!;
my %cutoff;
foreach my $speciesID(keys %sort){
	foreach my $rate(sort {$b <=> $a} keys %{$sort{$speciesID}}){
		if($rate >= $kmerCut){
			print NC "$speciesID\t".sprintf("%.2f\n", $rate*100);
			last;
		}else{
			if(!exists $cutoff{$speciesID} || $cutoff{$speciesID}{num} < $top){
				$cutoff{$speciesID}{num} += $sort{$speciesID}{$rate};
				$cutoff{$speciesID}{cut} = $rate;
			}else{
				last;
			}
		}
	}
}
close NC;

open IN,"$out" or die $!;
open OU,">$out.2" or die $!;
my %mark;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	# Extract the top N sequences with unique mapping confidence scores
	if(exists $cutoff{$seqid2speid{$t[2]}}){
		if($cutoff{$seqid2speid{$t[2]}}{cut} <= $t[-1]){
			if(!exists $mark{$seqid2speid{$t[2]}}){
				$mark{$seqid2speid{$t[2]}} = 1;
			}
			if($mark{$seqid2speid{$t[2]}} <= $top){
				print OU "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\n";
			}
			$mark{$seqid2speid{$t[2]}} += 1;
		}
	}elsif(exists $cutoff{$seqtid2S1taxid{$t[2]}}){
		if($cutoff{$seqtid2S1taxid{$t[2]}}{cut} <= $t[-1]){
			if(!exists $mark{$seqtid2S1taxid{$t[2]}}){
				$mark{$seqtid2S1taxid{$t[2]}} = 1;
			}
			if($mark{$seqtid2S1taxid{$t[2]}} <= $top){
				print OU "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\n";
			}
			$mark{$seqtid2S1taxid{$t[2]}} += 1;
		}
	}
}
close IN;
close OU;
`mv $out.2 $out`;

sub help{
print "
		Usage: perl $0

		--kout        kout file
		--top         [5]
		--db          Database version of PIP.
		--out         [./out.txt]
		--kmerCut     kmer cut off. [0.6875]

";
exit(0);
}

