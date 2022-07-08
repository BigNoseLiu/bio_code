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
GetOptions(
	"help"=>\$help,
	"af"=>\$freq_filter,
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
	#filter mutect auto filter flag
	if(	$arr[6] =~ /PASS/i ||
		$arr[6] =~ /^str_contraction$/i ||
		$arr[6] =~ /^clustered_events$/i ||
		$arr[6] =~ /^multiallelic$/i
	){
		my @tumor_formats = split(/:/,$arr[9]);
		my @tumor_afs = split(/,/,$tumor_formats[2]);
		my @tumor_genotype = split(/\//,$tumor_formats[0]);
		my @normal_formats = split(/:/,$arr[10]);
		my @normal_afs = split(/,/,$normal_formats[2]);
		#log low diff afs between tumor & normal
		my $is_lowDiffAF = 0;
		for( my$i=0;$i<@tumor_genotype;$i++ ){
			
			my $cmp_af = $tumor_afs[$i-1]/$normal_afs[$i-1];
			if( $tumor_genotype[$i]=~/1/i && $cmp_af < $min_odd ){
				$is_lowDiffAF = 1;
				print STDERR join("\t",$tumor_genotype[$i],$cmp_af,$i,$tumor_afs[$i-1],$normal_afs[$i-1])."\t";
			}
		}
		if($is_lowDiffAF == 1){
			$arr[6] .= ";liummLowDiffAF";
			print STDERR join("\t",@arr)."\n";
			#filter out low diff afs
			if( defined($freq_filter) ){
				next;
			}
		}
	}
	else{
		next;
	}
	if( $arr[0] eq $last_pos[0] && ($last_pos[1]-$arr[1]) < $max_dis && ($last_pos[1]-$arr[1]) > (0-$max_dis) ){
			$arr[6] .= ";liummLowDistance";
	}
	@last_pos = ($arr[0],$arr[1]);
	print join("\t",@arr)."\n";
}
