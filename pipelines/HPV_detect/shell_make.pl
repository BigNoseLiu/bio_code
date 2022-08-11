#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use FindBin qw($Bin $Script);

my %h_fqs = ();
my $line;
while( $line = <STDIN> ){
	#202205231446_210901006_B_POOLINHG18_ZMCAFY4085_L01_R2.fq.gz
	chomp $line;
	if( $line =~ /(POOLINHG\d+_[^_]+)_[^_]+_R1/ ){
		my $sample=$1;
		if( !defined($h_fqs{ $sample }{'R1'}) ){
			$h_fqs{ $sample }{'R1'} = $line;
		}
		else{
			print STDERR "dup sample for $sample R2\n";
		}
	}
	elsif( $line =~ /(POOLINHG\d+_[^_]+)_[^_]+_R2/ ){
		my $sample=$1;
		if( !defined($h_fqs{ $sample }{'R2'}) ){
			$h_fqs{ $sample }{'R2'} = $line;
		}
		else{
			print STDERR "dup sample for $sample R2\n";
		}
	}
	else{
		print STDERR "illegal fomrat for $line\n";
	}
}

my $out_dir = $ARGV[0];
`mkdir -p $out_dir`;
foreach my $sample( sort {$a cmp $b} keys(%h_fqs) ){
	print "sh $Bin/HPV_16_18.v1.sh $sample ".$h_fqs{$sample}{'R1'}." ".$h_fqs{$sample}{'R2'}." $out_dir/$sample\n";
}
