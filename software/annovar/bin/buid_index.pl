#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
my($in,$bins,$help);
GetOptions("i|in=s"=>\$in,"b|bins=i"=>\$bins,"help|?"=>\$help);
my$usage="";
#my$usage=<<<INFO;
#Usage:
 #       perl $0[options]
#Options:
#-i <file><in.input the database file for indexing>
#-b<INT><bin.size the size of the index's bin>
#INFO
die $usage if ($help || !$in || !$bins);
open IN,"$in" || die $!;
my $out = ${in} . ".idx";
my $size_in = -s $in;
my %idxLast=();
my %idxFirst=();
my $key_bin=0;
my $len = 0;
while(<IN>){
	if($_ =~ /^#/){
		$len+= length $_;
		next;
	}
	my($chr,$start,$end)=(split /\t/,$_)[0,1,2];
	$chr=~s/^chr//i;
	my$bin_start=$bins* int($start/$bins);
	$len+= length $_;
	$key_bin="$chr-$bin_start";
	if(exists $idxLast{$key_bin}){
		$idxLast{$key_bin}=$len;
	}
	else{
		$idxFirst{$key_bin}=$len- length $_;
		$idxLast{$key_bin}=$len;
	}
}
close IN;
open OUT,"> $out"||die$!;
print OUT "#BIN\t$bins\t$size_in\n";
foreach(sort keys %idxFirst){
	my($chr,$bin_start)= split /\-/,$_;
	next if($chr=~/^#/);
	print OUT "$chr\t$bin_start\t$idxFirst{$_}\t$idxLast{$_}\n";
}
close OUT;
