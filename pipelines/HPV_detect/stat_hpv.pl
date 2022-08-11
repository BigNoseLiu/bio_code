#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my $line;

my %h_sample = ();
open IN,"$Bin/hpv_sample_info.xls" or die $!;
while( $line = <IN> ){
	next if($line =~ /^#/);
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_sample{$arr[0]}{'mut'}=$arr[1];
}
close IN;

while(my $file_path = <STDIN> ){
	my $sample_id1 = "-";
	my $sample_id2 = "-";
	if( $file_path =~ /\/([^\.\/]+)_([^\._]+)\.[^\/]+$/ ){
		$sample_id1=$1."_".$2;
		$sample_id2=$2;
	}
	else{
		print STDERR "Err: Illeage file name for $file_path\n";
		exit;
	}

	open IN,"cat $file_path|" or die $!;
	my $unclass = 0;
	my $all_map = 0;
	my $hpv16_map = 0;
	my $hpv18_map = 0;
	while( $line = <IN> ){
		chomp $line;
		my @arr = split(/\t/,$line);
		if( $line =~ /^\s*\S+\s+(\S+)\s/ ){
			$arr[1] = $1;
		}
		$unclass = $arr[1] if($line =~ /unclassified/);
		$all_map = $arr[1] if($line =~ /root/);
		$hpv16_map = $arr[1] if($line =~ /Human papillomavirus type 16/);
		$hpv18_map = $arr[1] if($line =~ /Human papillomavirus type 18/);
		$h_sample{$sample_id2}{"seq"}{$sample_id1}{"total"} = $unclass+$all_map;
		$h_sample{$sample_id2}{"seq"}{$sample_id1}{"hpv16"} = $hpv16_map;
		$h_sample{$sample_id2}{"seq"}{$sample_id1}{"hpv18"} = $hpv18_map;
	}
	close IN;
}
print "#sample_id\tseq_id\thpv_type\tpredicted\tall_reads\thpv16_reads\thpv18_reads\n";
foreach my $sample_id2(sort {$a cmp $b} keys(%h_sample) ){
	if( !defined($h_sample{$sample_id2}{'mut'}) ){
		$h_sample{$sample_id2}{'mut'} = "-";
	}
	if( defined($h_sample{$sample_id2}{'seq'}) ){
		foreach my $sample_id1( sort {$a cmp $b} keys(%{$h_sample{$sample_id2}{"seq"}}) ){
			my $predict_result = "阴性";
			my $cut_off = 1000;
			if( $h_sample{$sample_id2}{"seq"}{$sample_id1}{"hpv16"} > $cut_off ){
				$predict_result = "HPV16阳性";
				if( $h_sample{$sample_id2}{"seq"}{$sample_id1}{"hpv18"} > $cut_off ){
					$predict_result .= "|HPV18阳性";
				}
			}
			elsif( $h_sample{$sample_id2}{"seq"}{$sample_id1}{"hpv18"} >= $cut_off ){
				$predict_result = "HPV18阳性";
			}
			print $sample_id2."\t$sample_id1\t".$h_sample{$sample_id2}{'mut'}."\t$predict_result\t".$h_sample{$sample_id2}{"seq"}{$sample_id1}{"total"}."\t".$h_sample{$sample_id2}{"seq"}{$sample_id1}{"hpv16"}."\t".$h_sample{$sample_id2}{"seq"}{$sample_id1}{"hpv18"}."\n";
		}
	}
	else{
		print $sample_id2."\t-\t".$h_sample{$sample_id2}{'mut'}."\n";
	}
}
