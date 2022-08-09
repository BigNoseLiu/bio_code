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
        Used to split vcf multiple alleles
Usage:  
        perl $0 [options] -in gatk.filter.vcf -ref hg19.fa -out gatk.filter.fix.vcf
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

my ($help,$input_file,$ref_fa_file,$outfile);
GetOptions(
	"help"=>\$help,
);
if (defined $help ) {
	&usageWithColor();
	exit 0;
}

my $line;

my %h_input = ();
my @arr_exomiser_yml = <STDIN>;
foreach $line( @arr_exomiser_yml ){
	foreach my $input( @ARGV ){
		next if( $input !~ /:/ );
		my @arr = split(/:/,$input);
		my $be_replace = shift @arr;
		my $to_replace = join( ":", @arr );
		$h_input{$be_replace} = $to_replace;
	}
}


#读取已录入的ped文件
my $fam = "FAM1";
my %h_ped1 = ();
if( defined($h_input{'liumm_ped_file_in'}) ){
	open IN,$h_input{'liumm_ped_file_in'} or die $!;
	while( $line = <IN> ){
		chomp $line;
		my @arr = split(/\t/,$line);
		$fam = $arr[0];
		$h_ped1{$arr[1]} = $line;
	}
	close IN;
}

my $all_merge_sample_id = "da_all_merge";
my %h_samples = ();
	open IN,$h_input{'liumm_vcf_file_in'} or die $!;
	open OUT,">".$h_input{'liumm_vcf_file_path'} or die $!;
	my @header_arr = ();
	while( $line = <IN> ){
		chomp $line;
		#vcf headers
		if( $line =~ /^#/ ){
			if( $line =~ /^#CHROM/ ){
				@header_arr = split(/\t/,$line);
				for( my$m=9;$m<scalar@header_arr;$m++ ){
					$header_arr[$m] =~ s/_dadup\d+$//;
					my $sample_id = $header_arr[$m];
					my $real_sample = $sample_id;
					
					$h_samples{$sample_id}++;
				}
				print OUT join("\t",@header_arr[0..8])."\t".join("\t",sort {$a cmp $b} keys(%h_samples))."\t$all_merge_sample_id\n";
			}
			else{
				print OUT $line."\n";
			}
			next;
		}
		my @arr = split(/\t/,$line);
		my %h_value = ();
		my $valule_all_merged = ".";
		for( my$m=9;$m<scalar@header_arr;$m++ ){
			#获取第一个非空的值，作为merge结果
			$valule_all_merged = $arr[$m] if($valule_all_merged =~ /^\./);
			if( defined($h_value{$header_arr[$m]}) && $h_value{$header_arr[$m]} =~ /^\d/ ){
			}
			else{
				$h_value{$header_arr[$m]} = $arr[$m];
			}
		}
		print OUT join("\t",@arr[0..8]);
		foreach my $sample_id(sort {$a cmp $b} keys(%h_samples)){
			print OUT "\t".$h_value{$sample_id};
		}
		print OUT "\t$valule_all_merged\n";
	}
	close IN;
	close OUT;


#修复ped文件，将没有的样本id加入
open OUT,">".$h_input{'liumm_ped_file_path'} or die $!;
foreach my $sample_id( sort {$a cmp $b} keys(%h_samples) ){
	if( defined($h_ped1{$sample_id}) ){
		print OUT $h_ped1{$sample_id}."\n";
	}
	else{
		print OUT "$fam\t$sample_id\t0\t0\tother\t0\n";
	}
}
print OUT "$fam\t$all_merge_sample_id\t0\t0\tother\t2\n";
close OUT;

#如果定义的proband不存在，取merge的结果
if( defined($h_input{'liumm_proband'}) && defined($h_samples{$h_input{'liumm_proband'}}) ){
}
else{
	$h_input{'liumm_proband'} = $all_merge_sample_id;
}

foreach $line( @arr_exomiser_yml ){
	foreach my $be_replace( 'liumm_proband','liumm_ped_file_path','liumm_vcf_file_path','liumm_hpo_ids','liumm_output_prefix'){
		my $to_replace = $h_input{$be_replace};
		$line =~ s/$be_replace/$to_replace/;
	}
	print $line;
}
