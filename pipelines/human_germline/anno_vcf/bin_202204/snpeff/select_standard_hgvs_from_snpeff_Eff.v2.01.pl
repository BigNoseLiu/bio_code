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

my %h_stand_trans = ();
open IN,"$Bin/transcript_info/GRCh37.p13_genomic.trans_info.txt" or die $!;
my $prior_c = 0;
while( my $line = <IN> ){
	chomp $line;
	$prior_c++;
	my @arr = split(/\t/,$line);
	for( my $i=1;$i<scalar@arr;$i+=2 ){
		if(!defined( $h_stand_trans{$arr[$i]})){
			$h_stand_trans{$arr[$i]}{'len'} = $arr[$i+1];
			$h_stand_trans{$arr[$i]}{'prior'} = $prior_c;
		}
	}
}
close IN;
#get snpeff annotated vcf file
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


my @lines = <STDIN>;
foreach my $line( @lines ){
	if( $line =~ /^#/ ){
		print $line;
		next;
	}
	#chomp $line;
	#$line = $line."\t";
	my $select_snpeff_anno = "-";
	if( $line =~ /ANN=([^;\t]+)[;\t]/ ){
		my $t_snpeff_annos = $1;
		my %h_impact = ();
		my %h_all_annos = ();
		foreach my $t_snpeff_anno( split(/,/, $t_snpeff_annos) ){
			my @t_arr = split(/\|/,$t_snpeff_anno);
			my $t_count = scalar@t_arr;
			for(my $i=$t_count;$i<16;$i++){
				$t_arr[$i] = "";
			}
			$h_impact{$t_arr[2]}{$t_arr[6]}{'prior'} =  100000000000000000;
			my $trans = $t_arr[6];
			$trans =~ s/\.\d+$//;
			if( defined( $h_stand_trans{$trans}{'prior'} ) ){
				$h_impact{$t_arr[2]}{$t_arr[6]}{'prior'} = $h_stand_trans{$trans}{'prior'};
				if( $t_arr[11]  eq "" ){
					$t_arr[11] = "-1/".$h_stand_trans{$trans}{'len'};
					if( $trans eq "NM_001146310" ){
						print STDERR "dfdfd\n";
					}
				}
			}
			$h_impact{$t_arr[2]}{$t_arr[6]}{'anno'} =  join('|',@t_arr);
		}
		my $t_snpeff_annos_new = "";
		foreach my $t_impact( 'HIGH','MODERATE','LOW','MODIFIER','' ){
			if( defined($h_impact{$t_impact}) ){
				foreach my $trans( sort {$h_impact{$t_impact}{$a}{'prior'} <=> $h_impact{$t_impact}{$b}{'prior'} } keys(%{$h_impact{$t_impact}}) ){
					$select_snpeff_anno = $h_impact{$t_impact}{$trans}{'anno'} if($select_snpeff_anno eq "-");
					$t_snpeff_annos_new .= ",".$h_impact{$t_impact}{$trans}{'anno'};
				}
			}
		}
		$t_snpeff_annos_new =~ s/^,//;
	print STDERR $line."\n";
		#$line =~ s/ANN=/abc/;
		#$line =~ s/ANN=$t_snpeff_annos/ANN_STANDARD=$select_snpeff_anno;ANN=$t_snpeff_annos_new/;
		$line =~ s/ANN=[^;\t]+([;\t])/ANN_STANDARD=$select_snpeff_anno;ANN=$t_snpeff_annos_new$1/;
	print STDERR $line."\n";
		exit;
		#$h_stand_trans{$arr[$i]}{'prior'} = $prior_c if(!defined( $h_stand_trans{$arr[$i]}{'prior'} ));
        	#2. effect最强的
		#my @t_arr_impact = sort sortImpact keys(%{$h_impact{$t_arr_trans_types[0]}});
        	#3. 转录本最长的
		#my @t_arr_trans = sort { $h_cdna_len{$b} <=> $h_cdna_len{$a} } keys(%{$h_impact{$t_arr_trans_types[0]}{$t_arr_impact[0]}});
		#$select_snpeff_anno = $h_impact{$t_arr_trans_types[0]}{$t_arr_impact[0]}{$t_arr_trans[0]};
	}
	#选好的注释，插入vcf Info
	#print $line;
}
