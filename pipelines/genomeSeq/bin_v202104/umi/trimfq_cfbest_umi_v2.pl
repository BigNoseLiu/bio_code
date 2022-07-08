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

my $fq_zip1 = $ARGV[0];
my $fq_zip2 = $ARGV[1];
my $out_fq1 = $ARGV[2];
my $out_fq2 = $ARGV[3];
my $out_stat1 = $ARGV[4];
my %h_stat = ();
my $line;
open IN,"$Bin/umi_code.cfbest.v2.txt" or die $!;
while( $line = <IN> ){ 
	chomp $line;
	$line =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
	$h_stat{$2}{"name"} = $3;
	$h_stat{$2}{"keep"} = $1;
	$h_stat{$2}{"all"} = 0;
	
}
close IN;
open IN1,"gunzip -c $fq_zip1|" or die $!;
open IN2,"gunzip -c $fq_zip2|" or die $!;
open OUT1,">$out_fq1" or die $!;
open OUT2,">$out_fq2" or die $!;
my $all_count = 0;
my $filter_count = 0;
while( $line = <IN1> ){
	my $out_str1 = $line;
	$line = <IN1>;
	$out_str1 .= $line;
	$line=~/^(.{7})/;
	my $umi1 = $1;	
	$line = <IN1>;
	$out_str1 .= $line;
	$line = <IN1>;
	$out_str1 .= $line;

	my $out_str2 = <IN2>;
	$line = <IN2>;
	$out_str2 .= $line;
	$line=~/^(.{7})/;
	my $umi2 = $1;	
	$line = <IN2>;
	$out_str2 .= $line;
	$line = <IN2>;
	$out_str2 .= $line;



	$h_stat{$umi1}{"all"}++;
	$h_stat{$umi2}{"all"}++;
	$h_stat{$umi1}{"fq1"}++;
	$h_stat{$umi2}{"fq2"}++;
	$all_count++;
	if( !defined($h_stat{$umi1}{"keep"}) || !defined($h_stat{$umi2}{"keep"}) || $h_stat{$umi1}{"keep"} !=1 || $h_stat{$umi2}{"keep"} != 1 ){
		next;
	}
	$filter_count++;
	$h_stat{$umi1}{"filter_fq1"}++;
	$h_stat{$umi2}{"filter_fq2"}++;
	print OUT1 $out_str1;
	print OUT2 $out_str2;
}
close IN1;
close IN2;
close OUT1;
close OUT2;

open OUT1, ">$out_stat1" or die $!;
print OUT1 "#all_count\t$all_count\n";
print OUT1 "#filter_count\t$filter_count\n";
my $filter_perc = $filter_count/$all_count*100;
$filter_perc =~ s/(\.\d\d).*$/$1/;
print OUT1 "#filter_perc(%)\t$filter_perc\n";
my @titles = ("fq1","fq2","filter_fq1","filter_fq2");
print OUT1 "#code\tname\t".join("\t",@titles)."\n";
foreach my $code( sort {$h_stat{$b}{"all"} <=> $h_stat{$a}{"all"}} keys(%h_stat) ){
	$h_stat{$code}{"name"} = "-" if(!defined($h_stat{$code}{"name"})  );
	$h_stat{$code}{"keep"} = 0 if(!defined($h_stat{$code}{"keep"})  );
	foreach my $title( @titles ){
		$h_stat{$code}{$title} = 0 if( !defined($h_stat{$code}{$title}) );
	}
	print OUT1 $code."\t".$h_stat{$code}{"keep"}."\t".$h_stat{$code}{"name"}."\t".$h_stat{$code}{"fq1"}."\t".$h_stat{$code}{"fq2"}."\t".$h_stat{$code}{"filter_fq1"}."\t".$h_stat{$code}{"filter_fq2"}."\n";
}
close OUT1;
