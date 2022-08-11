#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);



my $line;
my %h_accurate_mut = ();
my %h_mut = ();
my %h_mut_type = ();
open IN,"$Bin/mut_info.xls" or die $!;
while( $line = <IN> ){
	next if($line =~ /^#/);
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_mut_type{$arr[1]} = $arr[0];
	next if($arr[0] ne "SNV" && $arr[0] ne "indel" );
	$arr[5] = "-" if(!defined($arr[5]));
	$arr[6] = "-" if(!defined($arr[6]));
	$arr[7] = "-" if(!defined($arr[7]));
	$h_accurate_mut{$arr[2]}{join(":",$arr[5],$arr[6],$arr[7])}{$arr[1]}++;
	for( my $pos=$arr[3];$pos<=$arr[4];$pos++ ){
		$h_mut{$arr[2]}{$pos}{$arr[1]}++;
	}
}
close IN;

my %h_sample = ();
open IN,"$Bin/sample_info.xls" or die $!;
while( $line = <IN> ){
	next if($line =~ /^#/);
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_sample{$arr[0]}{'mut'}=$arr[2];
	if( defined($h_mut_type{$arr[2]}) ){
		$arr[2] = $h_mut_type{$arr[2]}."\t".$arr[2];
	}
	else{
		$arr[2] = "-\t".$arr[2];
	}
	$h_sample{$arr[0]}{'line'}=join("\t",@arr);
}
close IN;

print "#seq_id\tsample_id\tsample_type\tmut_type\tmut_name\tmut_info\taccurate_match\trough_match\taccurate_muts\taccurate_names\trough_muts\trough_names\n";
while(my $file_path = <STDIN> ){
	my $sample_id1 = "-";
	my $sample_id2 = "-";
	if( $file_path =~ /\/([^\.\/]+)_([^\._]+)\.[^\/]+$/ ){
		$sample_id1=$1."_".$2;
		$sample_id2=$2;
	}
	else{
		print STDERR "Err: Illeage file name for $file_path\n";
		exit;
	}

	my %h_accurate_stat = ();
	my %h_accurate_stat_mut = ();
	my %h_stat = ();
	my %h_stat_mut = ();
	open IN,"cat $file_path|" or die $!;
	while( $line = <IN> ){
		chomp $line;
		next if($line =~ /^#/);
		chomp $line;
		my @arr = split(/\t/,$line);
		my $t_mut = join(":",$arr[1],$arr[3],$arr[4]);
		if( defined($h_accurate_mut{$arr[0]}{$t_mut}) ){
			foreach my $mut(keys(%{$h_accurate_mut{$arr[0]}{$t_mut}})){
				my @arr_format = split(/:/,$arr[9]);
				$h_accurate_stat{ join(":",$arr[0],$arr[1],$arr[3],$arr[4],$arr_format[3],$arr_format[2]) }{$mut}++;
				$h_accurate_stat_mut{$mut}++;
			}
		}
		if( defined($h_mut{$arr[0]}{$arr[1]}) ){
			foreach my $mut(keys(%{$h_mut{$arr[0]}{$arr[1]}})){
				my @arr_format = split(/:/,$arr[9]);
				$h_stat{ join(":",$arr[0],$arr[1],$arr[3],$arr[4],$arr_format[3],$arr_format[2]) }{$mut}++;
				$h_stat_mut{$mut}++;
			}
		}
	}
	my $sample_info = "-\t-\t-\t-\t-";
	my $match = "-";
	my $accurate_match = "-";
	if( defined($h_sample{$sample_id2}) ){
		$sample_info = $h_sample{$sample_id2}{'line'};
		if( defined($h_accurate_stat_mut{ $h_sample{$sample_id2}{'mut'} }) ){
			$accurate_match = "accurate_yes";
		}
		else{
			$accurate_match = "accurate_no";
		}
		if( defined($h_stat_mut{ $h_sample{$sample_id2}{'mut'} }) ){
			$match = "rough_yes";
		}
		else{
			$match = "rough_no";
		}
	}
	#print $sample_id1."\t".$sample_info."\t".$accurate_match;
	my $accurate_all_vcf_info = "-";
	my $accurate_all_vcf_mut = "-";
	foreach my $vcf_info( sort {$a cmp $b} keys(%h_accurate_stat) ){
		my @muts = sort {$a cmp $b} keys(%{$h_accurate_stat{$vcf_info}});
		#print "\t$vcf_info\t".join("|",@muts);
		$accurate_all_vcf_info .= ";".$vcf_info;
		$accurate_all_vcf_mut .= ";".join("|",@muts);
	}
	#print "\t".$match;
	my $all_vcf_info = "-";
	my $all_vcf_mut = "-";
	foreach my $vcf_info( sort {$a cmp $b} keys(%h_stat) ){
		my @muts = sort {$a cmp $b} keys(%{$h_stat{$vcf_info}});
		#print "\t$vcf_info\t".join("|",@muts);
		$all_vcf_info .= ";".$vcf_info;
		$all_vcf_mut .= ";".join("|",@muts);
	}
	$accurate_all_vcf_info =~ s/^-;//;
	$accurate_all_vcf_mut =~ s/^-;//;
	$all_vcf_info =~ s/^-;//;
	$all_vcf_mut =~ s/^-;//;
	print $sample_id1."\t".$sample_info."\t$accurate_match\t$match\t$accurate_all_vcf_info\t$accurate_all_vcf_mut\t$all_vcf_info\t$all_vcf_mut\n";
	close IN;
}
