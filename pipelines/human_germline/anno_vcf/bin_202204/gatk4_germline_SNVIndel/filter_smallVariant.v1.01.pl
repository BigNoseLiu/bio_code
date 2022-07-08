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

my $max_pop = 0.05;
my @pop_freqs = ("gnomAD_genome_POPMAX","gnomAD_exome_POPMAX");
my $line;
#过滤hpo
my @hpo_ids = @ARGV;
my %h_result = ();
while( $line = <STDIN> ){
	chomp $line;
	if( $line =~ /^#/ ){
		print $line."\n";
		next;
	}
	my $pop_filter = 0;
	foreach my $info_title( @pop_freqs){
		if( $line =~ /$info_title=([^;]+);/ ){
			my $pop_freq = $1;
			if( $pop_freq ne '.' && $pop_freq > $max_pop ){
				$pop_filter = 1;
				#print STDERR $info_title."\t".$pop_freq."\n";
			}
		}
	}
	my $hpo_filter = 1;
	if( scalar@hpo_ids > 0 ){
		foreach my $hpo_id( @hpo_ids ){
			my $t_flag = 1;
			foreach my $t_hpo( split(/,/,$hpo_id) ){
				if( $line !~ /$t_hpo/ ){
					$t_flag = 0;
				}
			}
			$hpo_filter = 0 if($t_flag == 1);
		}
	}
	else{
		$hpo_filter = 0;
	}
	my @arr = split(/\t/,$line);
	if( $pop_filter == 0 && $hpo_filter == 0){
		if( $arr[9] !~ /\S/ ){
			$arr[9] = -1;
		}
		$h_result{$arr[9]}{$arr[0]."\t".$arr[1]."\t".$arr[2]."\t".$arr[3]}{$line}++;
	}

}
foreach my $score(sort {$b <=> $a} keys(%h_result)){
	foreach my $pos(sort {$a cmp $b} keys(%{$h_result{$score}})){
		foreach my $line(sort {$a cmp $b} keys(%{$h_result{$score}{$pos}})){
			print $line."\n";
		}
	}
}
