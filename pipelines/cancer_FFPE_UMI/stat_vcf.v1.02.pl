#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);



my $line;
my %h_all_muts = ();
my %h_accurate_mut = ();
my %h_mut_name = ();
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
	$h_mut_name{$arr[2]}{join(":",$arr[5],$arr[6],$arr[7])} = $arr[1];
	if(defined($arr[8]) && $arr[8] =~ /\S/){
		$h_mut_name{$arr[2]}{join(":",$arr[5],$arr[6],$arr[7])} = $arr[1]."($arr[8])";
	}
	$h_all_muts{ $h_mut_name{$arr[2]}{join(":",$arr[5],$arr[6],$arr[7])} } = 0;
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

my %h_stat_all_sample2muts = ();
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
	my $sample_info = $sample_id1."\t"."-\t-\t-\t-\t-";
	if( defined($h_sample{$sample_id2}) ){
		$sample_info = $sample_id1."\t".$h_sample{$sample_id2}{'line'};
	}

	$h_stat_all_sample2muts{$sample_info}{"match"} = 'no';
	$h_stat_all_sample2muts{$sample_info}{"match"} = '-' if( defined($h_sample{$sample_id2}) && defined($h_sample{$sample_id2}{'mut'}) && $h_sample{$sample_id2}{'mut'} eq '阴性' );

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
		my $aa_change="-";
		if( $line =~ /Func.refGeneWithVer=([^;]+);/ ){
			$aa_change=$1;
		}
		if( $line =~ /AAChange.refGeneWithVer=([^;]+);/ ){
			my $t_str = $1;
			$aa_change=$t_str if( $t_str ne '.' );
		}
		my @arr_format = split(/:/,$arr[9]);
		if( defined($h_accurate_mut{$arr[0]}{$t_mut}) ){
			$h_stat_all_sample2muts{$sample_info}{"var_freq"}{ $h_mut_name{$arr[0]}{$t_mut} } = $arr_format[2]*100;
			$h_all_muts{ $h_mut_name{$arr[0]}{$t_mut} }++;
			foreach my $mut(keys(%{$h_accurate_mut{$arr[0]}{$t_mut}})){
				$h_accurate_stat{ join(":",$aa_change,$arr[0],$arr[1],$arr[3],$arr[4],$arr_format[3],$arr_format[2]) }{$mut}++;
				$h_accurate_stat_mut{$mut}++;
			}
		}
		#if( defined($h_mut{$arr[0]}{$arr[1]}) || ( defined( $h_mut_name{$arr[0]}{$t_mut} ) && $h_mut_name{$arr[0]}{$t_mut} =~ /MET/ && $line =~ /exon14/ )){
		if( defined($h_mut{$arr[0]}{$arr[1]}) ){
			foreach my $mut(keys(%{$h_mut{$arr[0]}{$arr[1]}})){
				$h_stat{ join(":",$aa_change,$arr[0],$arr[1],$arr[3],$arr[4],$arr_format[3],$arr_format[2]) }{$mut}++;
				$h_stat_mut{$mut}++;
			}
		}
	}
	close IN;

	my $match = "-";
	my $accurate_match = "-";
	if( defined($h_sample{$sample_id2}) ){
		if( defined($h_accurate_stat_mut{ $h_sample{$sample_id2}{'mut'} }) ){
			$accurate_match = "accurate_yes";
			$h_stat_all_sample2muts{$sample_info}{"match"} = 'yes';
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
	print $sample_info."\t$accurate_match\t$match\t$accurate_all_vcf_info\t$accurate_all_vcf_mut\t$all_vcf_info\t$all_vcf_mut\n";
}
open OUT,">".$ARGV[0] or die $!;
my @arr_all_muts = sort {$a cmp $b } keys(%h_all_muts);
#my @arr_all_muts = sort {$h_all_muts{$b} <=> $h_all_muts{$a} } keys(%h_all_muts);
print OUT "#seq_id\tsample_id\tsample_type\tmut_type\tmut_name\tmut_info\tmatch\t".join("\t",@arr_all_muts )."\n";
foreach my $sample(sort {$a cmp $b} keys(%h_stat_all_sample2muts)){
	print OUT $sample."\t".$h_stat_all_sample2muts{$sample}{"match"};
	foreach my $mut( @arr_all_muts ){
		print OUT "\t";
		if( defined($h_stat_all_sample2muts{$sample}{"var_freq"}{$mut}) ){
			print OUT $h_stat_all_sample2muts{$sample}{"var_freq"}{$mut};
		}
	}
	print OUT "\n";
}
close OUT
