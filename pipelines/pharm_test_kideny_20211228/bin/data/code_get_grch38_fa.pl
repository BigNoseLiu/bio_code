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


my $line;
my %h_target = ();
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_target{$arr[4]}{$arr[5]}{$arr[3]}++;
	#print $arr[4]."\t".$arr[5]."\t".$arr[0]."\n";
}



my $ref_files = "../../../../database/ref_genomes/GRCh38.d1.vd1/GRCh38.d1.vd1.fa";

foreach my $t_ref( split(/,/,$ref_files) ){
	open REF,"$t_ref" or die "Fail opening $t_ref\n";
	my $line;
	$/ = '>';
	my $chr = "";
	while( $line = <REF> ){
		chomp $line;
		next if( $line =~ /^\s*$/ );
		if( $line =~ s/^(\S+)\s[^\n]*\n// ){
			$chr = $1;
		}
		else{
			print STDERR "wrong formated fa:$t_ref\n";
			next;
		}
		last if($chr =~ /_/);
		$line =~ s/[>\s]//g;
		foreach my $pos(sort {$a <=> $b} keys(%{$h_target{$chr}})){
			print ">$chr:$pos:".join(",",sort {$a cmp $b} keys(%{$h_target{$chr}{$pos}}))."\n";
			print substr($line,$pos-199,400)."\n";
		}
	}
	close REF;
}
