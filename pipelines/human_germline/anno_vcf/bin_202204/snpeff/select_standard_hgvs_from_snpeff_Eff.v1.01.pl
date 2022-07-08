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
        Used to 选择snpeff注释的最合适的基因转录本hgvs【1. NM_或NR的（若没有，则随便）.  2. effect最强的   3. 转录本最长的】
	note：vcf必须是split multiallelic sites into biallelic records【done by 'bcftools norm'】
Usage:  
        perl $0 <snpeff_anno.vcf >snpeff_anno.select.vcf
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
#get snpeff annotated vcf file
my @lines = <STDIN>;
my %h_cdna_len = ();
foreach my $line( @lines ){
	if( $line =~ /^#/ ){
		if( $line =~ /^##INFO=<ID=ANN,/ ){
			print '##INFO=<ID=ANN_STANDARD,Number=.,Type=String,Description="representation for ANN anno">'."\n";
		}
		print $line;
		next;
	}
	if( $line =~ /ANN=([^;\t]+)[;\s]/ ){
		my $t_snpeff_annos = $1;
		foreach my $t_snpeff_anno( split(/,/, $t_snpeff_annos) ){
			my @t_arr = split(/\|/,$t_snpeff_anno);
			$t_arr[11] = "-" if( !defined($t_arr[11]) );
			$t_arr[6] = "-" if( !defined($t_arr[6]) );
			if( $t_arr[11] =~ /\/(\d+)$/ ){
				$h_cdna_len{$t_arr[6]} = $1;
			}
			elsif( !defined($h_cdna_len{$t_arr[6]}) ){
				$h_cdna_len{$t_arr[6]} = 1;
			}
		}

	}
}

sub sortImpact{
	my $a1 = 6;
	my $b1 = 6;
	$a1 = 1 if( $a =~ /^HIGH$/i );
	$b1 = 1 if( $b =~ /^HIGH$/i );
	$a1 = 2 if( $a =~ /^MODERATE$/i );
	$b1 = 2 if( $b =~ /^MODERATE$/i );
	$a1 = 3 if( $a =~ /^LOW$/i );
	$b1 = 3 if( $b =~ /^LOW$/i );
	$a1 = 3 if( $a =~ /^MODIFIER$/i );
	$b1 = 3 if( $b =~ /^MODIFIER$/i );
	return $a1 <=> $b1;
}

sub sortTransType{
	my $a1 = 3;
	my $b1 = 3;
	$a1 = 1 if( $a =~ /^nm$/i );
	$b1 = 1 if( $b =~ /^nm$/i );
	$a1 = 2 if( $a =~ /^nr$/i );
	$b1 = 2 if( $b =~ /^nr$/i );
	$a1 = 4 if( $a =~ /^xm$/i );
	$b1 = 4 if( $b =~ /^xm$/i );
	$a1 = 5 if( $a =~ /^xr$/i );
	$b1 = 5 if( $b =~ /^xr$/i );
	return $a1 <=> $b1;
}

foreach my $line( @lines ){
	if( $line =~ /^#/ ){
		next;
	}
	chomp $line;
	$line = $line."\n";
	my $select_snpeff_anno = "-";
	if( $line =~ /ANN=([^;\t]+)[;\s]/ ){
		my $t_snpeff_annos = $1;
		my %h_impact = ();
		foreach my $t_snpeff_anno( split(/,/, $t_snpeff_annos) ){
			my @t_arr = split(/\|/,$t_snpeff_anno);
			my $t_trans_type = "-";
			$t_arr[2] = "-" if( !defined($t_arr[2]) );
			$t_arr[6] = "-" if( !defined($t_arr[6]) );
			if( $t_arr[6] =~ /^(..)/ ){
				$t_trans_type = $1;
				$t_trans_type = "-" if( $t_trans_type !~ /nm/i && $t_trans_type !~ /nr/i && $t_trans_type !~ /xm/i && $t_trans_type !~ /xr/i );
			}
			$h_impact{$t_trans_type}{$t_arr[2]}{$t_arr[6]} =  $t_snpeff_anno ;
		}
        	#1. NM_或NR的（若没有，则随便）.
		my @t_arr_trans_types = sort sortTransType keys(%h_impact);
        	#2. effect最强的
		my @t_arr_impact = sort sortImpact keys(%{$h_impact{$t_arr_trans_types[0]}});
        	#3. 转录本最长的
		my @t_arr_trans = sort { $h_cdna_len{$b} <=> $h_cdna_len{$a} } keys(%{$h_impact{$t_arr_trans_types[0]}{$t_arr_impact[0]}});
		$select_snpeff_anno = $h_impact{$t_arr_trans_types[0]}{$t_arr_impact[0]}{$t_arr_trans[0]};
	}
	#选好的注释，插入vcf Info
	$line =~ s/ANN=/ANN_STANDARD=$select_snpeff_anno;ANN=/;
	print $line;
}
