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
        Used to fix snpeff impact name
Usage:  
        perl $0 <snpeff_anno.vcf >snpeff_anno.fix.vcf
        Options:
	     -help : help info
Author & Contact:
	Mingming Liu
Last updated:
        2019-1-10
USAGE
	print color "reset";#change back the text color
}

my ($help);
GetOptions(
	"help"=>\$help,
);
if (defined $help ) {
	&usageWithColor();
	exit 0;
}

sub err_print{
	print STDERR color "red";#change the text color
	foreach my $t_cmd( @_ ){
		print STDERR $t_cmd."\n";
	}
	print STDERR color "reset";#change back the text color
}
my $line;
while( $line = <STDIN> ){
	if( $line =~ /^#/ ){
		print $line;
		next;
	}
	if( $line =~ /ANN=([^;]+)[;\s]/ ){
		my $t_snpeff_annos = $1;
		my $fix_snpeff_annos = "-";
		foreach my $t_snpeff_anno( split(/,/, $t_snpeff_annos) ){
			my @t_arr = split(/\|/,$t_snpeff_anno);
			#fix the effect type
			if( $t_arr[1] =~ /exon_loss_variant/ ){
				#1. exon_loss_variant
				$t_arr[1] = "exon_loss_variant";
			}
			elsif( $t_arr[1] =~ /frameshift_variant/i ){
				#2. frameshift_variant
				$t_arr[1] = "frameshift_variant";
			}
			elsif( $t_arr[1] =~ /splice_donor_variant/i || $t_arr[1] =~ /splice_acceptor_variant/i ){
				#3. splice_donor_variant or splice_acceptor_variant
				$t_arr[1] =~ s/intron_variant&//i;	$t_arr[1] =~ s/&intron_variant//i;
				$t_arr[1] =~ s/splice_region_variant&//i;	$t_arr[1] =~ s/&splice_region_variant//i;
			}
			#last. sort eff
			my @t_eff_arr = split(/&/,$t_arr[1]);
			$t_arr[1] = join("&", ( sort {$a cmp $b} @t_eff_arr ) );

			$fix_snpeff_annos .= ",".join('|',@t_arr);
		}
		$fix_snpeff_annos =~ s/^-,//;
		$line =~ s/$t_snpeff_annos/$fix_snpeff_annos/;
	}
	print $line;
}
