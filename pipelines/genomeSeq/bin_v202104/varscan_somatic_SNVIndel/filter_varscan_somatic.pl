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

sub usageWithColor {
	print color "green";#change the text color
	print <<USAGE;
Description:
        Used to ...
Usage:  
        perl $0 [options] <gatk.auto_filter.input.vcf     >gatk.self_filter.output.vcf
        Options:
	     -help : reveal help info
	     -af : filter low af diff alts if defined
	     -out  <str>  : output file path
        Example:
             perl $0 -out  output_file
Author & Contact:
	Mingming Liu
Last updated:
        2013-10-16
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}

my ($help,$freq_filter,$outfile);
my $min_coverage = 5;
GetOptions(
	"help"=>\$help,
	"af"=>\$freq_filter,
	"min-cov=i"=>\$min_coverage,
	"out=s"=>\$outfile,
);
if (defined $help ) {
	&usageWithColor();
	exit 0;
}

#starting
my $line;
my @last_pos = ("-",-1000);
my $min_odd = 1.25;
my $max_dis = 10;
while( $line = <STDIN> ){
	chomp $line;

	#vcf headers
	if( $line =~ /^#/ ){
		print $line."\n";
		next;
	}

	my @arr = split(/\t/,$line);
	$arr[7] .= ";";
	next if( $arr[7] !~ /somatic/i );
	if( $arr[7] =~ /SPV=([^;]+);/ ){
		my $p = $1;
		next if($p>0.05);
	}
	my @tumor_formats = split(/:/,$arr[10]);
	my @tumor_dp4 = split(/,/,$tumor_formats[6]);
	my @normal_formats = split(/:/,$arr[9]);
	if( ($tumor_dp4[2] + $tumor_dp4[3]) < $min_coverage ){
		next;
	}
	print join("\t",@arr)."\n";
}
