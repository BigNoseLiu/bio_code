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
        Used to 过滤cnvkit结果cns文件
Usage:  
        perl $0 [options] <cnvkit.cns  >cnvkit.filtered.cns
        Options:
	     -help : reveal help info
	     -odds  <str>  : 常染色体Max_odds,常染色体Min_odds,X染色体Max_odds,X染色体Min_odds,Y染色体Max_odds,Y染色体Min_odds
        Example:
             perl $0 -odds 0.3,-0.4,0.3,-0.4,0.3,-0.4 <B1701.markDup.sorted.cns  >B1701.markDup.sorted.filtered.cns
Author & Contact:
	Mingming Liu
Last updated:
        2018-11-16
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}

my ($help, $outfile );
my $limit_odds = "0.3,-0.4,0.8,-3,0.8,-3";
GetOptions(
	"help"=>\$help,
	"out=s"=>\$outfile,
	"odds=s"=>\$limit_odds,
);

if (defined $help ) {
	&usageWithColor();
	exit 0;
}

#获取常染色体，X染色体，Y染色体的log2比例阈值
my ($auto_max,$auto_min,$x_max,$x_min,$y_max,$y_min) = split(/,/,$limit_odds);

#header
my $line = <STDIN>;
print $line;
while( $line = <STDIN> ){
	chomp $line;
	if( $line =~ /^#/ || $line =~ /^\s*$/){
		print $line."\n";
		next;
	}
	my @arr = split(/\t/,$line);
	my ($max_odds, $min_odds) = ($auto_max, $auto_min);
	($max_odds, $min_odds) = ($x_max, $x_min) if( $arr[0] =~ /x/i );
	($max_odds, $min_odds) = ($y_max, $y_min) if( $arr[0] =~ /y/i );
	if( $arr[4] > $max_odds || $arr[4] < $min_odds ){
		print $line."\n";
	}
}
