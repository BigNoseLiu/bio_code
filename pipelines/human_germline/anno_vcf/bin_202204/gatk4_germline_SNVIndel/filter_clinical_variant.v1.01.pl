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
        Used to filter clinical variant
Usage:  
        perl $0 [options]  gatk.hg19.merge.txt >gatk.filter.vcf
        Options:
	     -help : reveal help info
        Example:
        	perl $0 <gatk.in.vcf  >gatk.out.vcf
Author & Contact:
	Mingming Liu
Last updated:
        2019-1-14
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}

my ($help,$ref_fa_file,$outfile);
GetOptions(
	"help"=>\$help,
);
if (defined $help ) {
	&usageWithColor();
	exit 0;
}

my $line;
#加载exomiser结果

#加载annovar.txt结果，获取纯合、杂合类型
my $annovar_txt = $ARGV[0] or die $!;
open IN, $annovar_txt or die $!;

#获取标题所对应的列序号
$line = <IN>;	chomp $line;
print $line."\n";
my @arr_header = split(/\t/,$line);
my %h_header = ();
for(my$i=0;$i<@arr_header;$i++){
	my $t_header = $arr_header[$i];
	$h_header{$t_header} = $i;
}

my $max_freq = 0.05;

while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	my $freq_filter = 0;
	foreach my $af_header( keys(%h_header) ){
		foreach my $db( 'ExAC','1000g2015aug','gnomadExome', 'gnomadGenome' ){
			if( $af_header =~ /^$db\_(\S+)$/ ){
				my $pop_type = $1;
				$arr[$h_header{$af_header}] = 0 if( $arr[$h_header{$af_header}] eq "." );
				if( $pop_type =~ /all/i || $pop_type =~ /eas/i ){
					$freq_filter = 1 if( $arr[$h_header{$af_header}] > $max_freq );
				}
			}
		}
	}
	if( $freq_filter == 1 ){
		next;
	}
	print $line."\n";
}

close IN;
