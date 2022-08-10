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
        Used to extract the reference sequence.
Reference genomes:
	/ifs1/NGB_DB/liumingming/XArchive4Gene/database/ftp_oriented/virus/All/virus.all.fa
	/ifs1/NGB_DB/liumingming/XArchive4Gene/database/ftp_oriented/genome/human/hg19/hg19.fa
Usage:  
        perl $0 [options]
        Options:
	     -help : reveal help info
	     -out  <str>  : output file path
        Example:
             perl $0 -whole -chr chr1 -pos 23223:200 -ref hg19.fa
             perl $0 -chr hr1 -pos 23223:200 -ref hg19.fa
             perl $0 -num -chr 1 -pos 23223:200 -ref hg19.fa
Author & Contact:
	Mingming Liu
        liumingming\@genomics.cn
Last updated:
        2014-2-27
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}

my ($help,$not_fuzzy_match,$number_chr,$ref_files,$chr,$pos);
#default settings
$ref_files = "/ifs1/NGB_DB/liumingming/XArchive4Gene/database/ftp_oriented/genome/human/hg19/hg19.fa";
GetOptions(
	"help"=>\$help,
	"whole"=>\$not_fuzzy_match,
	"number"=>\$number_chr,
	"ref=s"=>\$ref_files,
	"chr=s"=>\$chr,
	"pos=s"=>\$pos,
);
if (defined $help || (!defined $ref_files) || (!defined $chr) || (!defined $pos) ||
	(defined( $number_chr ) && $chr =~ /\D/ ) || $pos !~ /^\d+:\d+$/
) {
	&usageWithColor();
	exit 0;
}

foreach my $t_ref( split(/,/,$ref_files) ){
	open REF,"$t_ref" or die "Fail opening $t_ref\n";
	my $line;
	$/ = '>';
	my $title = "";
	my $chr_number = 0;
	while( $line = <REF> ){
		chomp $line;
		next if( $line =~ /^\s*$/ );
		if( $line =~ s/^([^\n]*)\n// ){
			$title = $1;
			$chr_number ++;
		}
		else{
			print STDERR "wrong formated fa:$t_ref\n";
			next;
		}
		$line =~ s/[>\s]//g;
		if( ( !defined( $number_chr ) && !defined( $not_fuzzy_match ) && $title !~ /$chr/) ||
		    ( defined( $number_chr ) && $chr_number != $chr ) ||
		    ( defined( $not_fuzzy_match ) && $title ne $chr )
		){
			next;
		}
		my ($start_pos,$len) = split( /:/,$pos );
		print ">$title\n";
		#print ">$title:".join(" ",@ARGV)."\n";
		print substr($line,$start_pos-1,$len)."\n";
	}
	close REF;
}
