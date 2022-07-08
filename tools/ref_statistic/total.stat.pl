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
        perl $0 [options]
        Options:
	     -help : reveal help info
	     -out  <str>  : output file path
        Example:
             perl $0 -out  output_file
Author & Contact:
	Mingming Liu
        liumingming\@genomics.cn
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

my ($help,$in_dir,$out);
GetOptions(
	"help"=>\$help,
	"in=s"=>\$in_dir,
	"out=s"=>\$out,
);
if (defined $help || (!defined $out) ) {
	&usageWithColor();
	exit 0;
}
my $line;
my %h_total;
my ($total,$n_base,$gc_base,$lc_base,$n_count) = ("total","n_base","gc_base","lc_base","n_count");
#chr    base    continous_len('-' stand for other meaning)      count
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split( /\t/,$line );
	if($arr[2] eq "-"){
		$h_total{ $arr[0] }{ $total } += $arr[3];
		if( $arr[1] =~ /[atcg/ ){
			$h_total{ $arr[0] }{ $lc_base } += $arr[3];
		}
		elsif( $arr[1] =~ /[cCgG]/ ){
			$h_total{ $arr[0] }{ $gc_base} += $arr[3];
		}
		elsif( $arr[1] =~ /[nN]/ ){
			$h_total{ $arr[0] }{ $n_base } += $arr[3];
		}
	}
}

open OUT,'>',"$out" or die "Fail opening $out\n";
print OUT "Start time:".time."\n";
close OUT;
