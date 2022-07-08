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

my ( $filter_flag );
GetOptions(
	"filter"=>\$filter_flag,
);

my $vcf1 = $ARGV[0];
my $vcf2 = $ARGV[1];
open IN1,"$vcf1" or die $!;
open IN2,"$vcf2" or die $!;
my $all_count = 0;
my $filter_count = 0;
my %h_vcf1 = ();
my $count1 = 0;
my $count2 = 0;
my $count_common = 0;
my $line;
while( $line = <IN1> ){
	next if( $line =~ /^#/ );
	next if( defined($filter_flag) && $line !~ /PASS/ );
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_vcf1{$arr[0]."\t",$arr[1]."\t".$arr[3]."\t".$arr[4]}++;
	$count1++;
}

while( $line = <IN2> ){
	next if( $line =~ /^#/ );
	next if( defined($filter_flag) && $line !~ /PASS/ );
	chomp $line;
	my @arr = split(/\t/,$line);
	if(defined($h_vcf1{$arr[0]."\t",$arr[1]."\t".$arr[3]."\t".$arr[4]})){
		$count_common++;
	}
	$count2++;
}

close IN1;
close IN2;
my $ratio1 = $count_common/$count1*100;
$ratio1 =~ s/(\.\d\d).*$/$1/;
my $ratio2 = $count_common/$count2*100;
$ratio2 =~ s/(\.\d\d).*$/$1/;
$vcf1 =~ s/.*\/([^\/]+)$/$1/;
$vcf2 =~ s/.*\/([^\/]+)$/$1/;
print "$count_common\t$ratio1\t$ratio2\t$count1\t$count2\t$vcf1\t$vcf2\n";
