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

my %h_stat = ();
my %h_type = ();
my $line;

while( $line = <STDIN> ){
	if( $line =~ /^#/ ){
		next;
	}
	chomp $line;
	$line = $line."\n";
	my $sig = "-";
	if( $line =~ /CLNSIG=([^;]+);/ ){
		$sig = $1;
	}
	if( $line =~ /ANN=([^;\t]+)[;\s]/ ){
		my $t_snpeff_annos = $1;
		my %h_impact = ();
		foreach my $t_snpeff_anno( split(/,/, $t_snpeff_annos) ){
			my @t_arr = split(/\|/,$t_snpeff_anno);
			$t_arr[2] = "-" if( !defined($t_arr[2]) );
			$t_arr[6] = "-" if( !defined($t_arr[6]) );
			foreach my $t_type(split(/&/,$t_arr[1])){
				$h_stat{$t_arr[3]}{$t_arr[6]}{$t_type}{$sig}++;
			}
		}
	}
}


foreach my $gene( sort {$a cmp $b} keys(%h_stat) ){
	foreach my $trans( sort {$a cmp $b} keys(%{$h_stat{$gene}}) ){
		foreach my $type( sort {$a cmp $b} keys(%{$h_stat{$gene}{$trans}}) ){
			foreach my $sig( sort {$a cmp $b} keys(%{$h_stat{$gene}{$trans}{$type}}) ){
				print "$gene\t$trans\t$type\t$sig\t".$h_stat{$gene}{$trans}{$type}{$sig}."\n";
			}
		}
	}
}
