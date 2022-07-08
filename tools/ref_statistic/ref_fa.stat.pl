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
        Used to statiistic the structure of reference genome.
Usage:  
        Options:
	     -help : reveal help info
	     -ref  <str>  : reference fasta file name
	     -out  <str>  : output file name
        Example:
             perl $0 -ref hg19.fa -out output_file
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

my ( $help, $out );
my ( $ref_fa_file );
GetOptions(
	"help"=>\$help,
	"out=s"=>\$out,
	"ref=s"=>\$ref_fa_file,
);
if( defined($help) || !defined($ref_fa_file) || !defined($out) ) {
	&usageWithColor();
	exit 0;
}
###################################
##description here               ##
###################################

#echo the start time
print "perl $0 -ref $ref_fa_file -out $out. From Start-time: ".localtime()." To ";

#no meaning str
my $NMS = "title_liumm";
my $line;
my (%h_sta,%h_base);
my $chr_name = "";
my $last_base = $NMS;
my $current_base = "";
my $continous_len = 0;

open REF,"$ref_fa_file" or die "Fail opening $ref_fa_file\n";
while( $line = <REF> ){
	chomp $line;
	if( $line =~ /^>(.*)$/ ){
		$chr_name = $1;
		$last_base = $NMS;
		next;
	}
	while( $line =~ /(.)/g ){
		$current_base = $1;
		$h_base{$chr_name}{$current_base}++;
		if( $last_base ne uc($current_base) ){
			$h_sta{$chr_name}{$last_base}{$continous_len}++;
			$last_base = uc($current_base);
			$continous_len = 1;
		}
		else{
			$continous_len++;
		}
	}
}
close REF;
$h_sta{$chr_name}{$last_base}{$continous_len}++;

open OUT,'>',"$out" or die "Fail writing $out\n";
#title
print OUT "#chr\tbase\tcontinous_len('-' stand for other meaning)\tcount\n";
#keys
foreach my $chr( sort {$a cmp $b} keys(%h_base) ){
	#2-d hash 
	foreach my $base( sort {$a cmp $b} keys( %{$h_base{$chr}} ) ){
		print OUT $chr."\t".$base."\t-\t".$h_base{$chr}{$base}."\n";
	}
	foreach my $base( sort {$a cmp $b} keys( %{$h_sta{$chr}} ) ){
		next if( $base eq $NMS );
		foreach my $len( sort {$a <=> $b} keys( %{$h_sta{$chr}{$base}} ) ){
			print OUT $chr."\t".$base."\t".$len."\t".$h_sta{$chr}{$base}{$len}."\n";
		}
	}
}
close OUT;

#echo the end time
print " End-time: ".localtime()."\n";
