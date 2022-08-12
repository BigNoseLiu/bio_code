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
        	perl $0 <vcf_file_path  >merge.out.vcf
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
my %h_vcf_header = ();
my @main_chr = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT');
my %h_chrs = ();
my %h_header_chrs = ();
my %h_result= ();
my %h_format= ();
my @main_format= ('GT','AD','AF','DP','FILT');
my %h_samples = ();
while(my $vcf_file=<STDIN>){
	chomp $vcf_file;
	open IN,$vcf_file or die $!;
	my @header_arr = ();
	while( $line = <IN> ){
		chomp $line;
		#vcf headers
		if( $line =~ /^#/ ){
			if( $line =~ /^##([^=]+)=/ ){
				my $header_type = $1;
				if( $line =~ /^##contig=<ID=([^,]+),/ ){
					my $t_chr = $1;
					if(!defined($h_chrs{$t_chr}{'info'})){
						$h_chrs{$t_chr}{'info'} = $line;
						$h_chrs{$t_chr}{'prior'} = 0;
					}
				}
				else{
					$h_vcf_header{$header_type}{$line}++;
				}
			}
			elsif( $line =~ /^#CHROM/ ){
				@header_arr = split(/\t/,$line);
				for( my$m=9;$m<scalar@header_arr;$m++ ){
					my $sample_id = $header_arr[$m];
					if( defined( $h_samples{$sample_id} ) ){
						$h_samples{$sample_id}++;
						$sample_id .= "_dadup".$h_samples{$sample_id};
					}
					$h_samples{$sample_id}++;
					$header_arr[$m] = $sample_id;
				}
			}
			next;
		}
		my @arr = split(/\t/,$line);
		$h_chrs{$arr[0]}{'prior'} = 0;
		my @format_arr = split(/:/,$arr[8]);
		$h_result{$arr[0]}{$arr[1]}{$arr[3]}{$arr[4]}{"info"} = join("\t",@arr[0..7]) if(!defined( $h_result{$arr[0]}{$arr[1]}{$arr[3]}{$arr[4]} ));
		for( my$m=9;$m<scalar@arr;$m++ ){
			#获取当前样本检测结果
			my @format_value_arr = split(/:/,$arr[$m]);
			for (my $i=0;$i<scalar@format_arr;$i++ ){
				$h_format{$format_arr[$i]} = 0;
				$h_result{$arr[0]}{$arr[1]}{$arr[3]}{$arr[4]}{"format"}{$header_arr[$m]}{$format_arr[$i]} = $format_value_arr[$i];
			}
			$h_result{$arr[0]}{$arr[1]}{$arr[3]}{$arr[4]}{"format"}{$header_arr[$m]}{'FILT'} = $arr[6];
			$h_format{'FILT'} = 0;
		}
	}
	close IN;
}

#chr排序
for (my $i=0;$i<scalar@main_chr;$i++ ){
	if( defined($h_chrs{$main_chr[$i]}{'prior'}) ){
		$h_chrs{$main_chr[$i]}{'prior'} = $i+1;
	}
}
my $count_chr = 1000;
foreach my $chr(sort {$a cmp $b} keys(%h_chrs)){
	if( $h_chrs{$chr}{'prior'} == 0  ){
		$h_chrs{$chr}{'prior'} = $count_chr;
		$count_chr++;
	}
}

#format字段排序
for (my $i=0;$i<scalar@main_format;$i++ ){
	if( defined($h_format{$main_format[$i]}) ){
		$h_format{$main_format[$i]} = $i+1;
	}
}
my $count_format= 1000;
foreach my $format(sort {$a cmp $b} keys(%h_format)){
	if( $h_format{$format} == 0  ){
		$h_format{$format} = $count_format;
		$count_format++;
	}
}
#打印headerss
$h_vcf_header{'FILT'}{'##FORMAT=<ID=FILT,Number=R,Type=Integer,Description="filter tag for all samples">'}++;
foreach my $header_type( sort {$a cmp $b} keys(%h_vcf_header)){
	foreach my $line( sort {$a cmp $b} keys(%{$h_vcf_header{$header_type}}) ){
		print $line."\n";
	}
	if( $header_type eq "INFO" ){
		foreach my $chr(sort {$h_chrs{$a}{'prior'} <=> $h_chrs{$b}{'prior'}} keys(%h_chrs)){
			if(defined($h_chrs{$chr}{'info'})){
				print $h_chrs{$chr}{'info'}."\n";
			}
		}
	}
}

my @arr_samples = sort {$a cmp $b} keys(%h_samples);
my @arr_formats = sort {$h_format{$a} <=> $h_format{$b}} keys(%h_format);
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	".join("\t",@arr_samples)."\n";
foreach my $chr(sort {$h_chrs{$a}{'prior'} <=> $h_chrs{$b}{'prior'} } keys(%h_chrs)){
	foreach my $pos( sort {$a <=> $b} keys(%{$h_result{$chr}}) ){
		foreach my $ref( sort {$a cmp $b} keys(%{$h_result{$chr}{$pos}}) ){
			foreach my $alt( sort {$a cmp $b} keys(%{$h_result{$chr}{$pos}{$ref}}) ){
				print $h_result{$chr}{$pos}{$ref}{$alt}{"info"}."\t".join(":",@arr_formats);
				foreach my $sample_id( @arr_samples ){
					my $sample_format = "";
					foreach my $format( @arr_formats ){
						if(defined($h_result{$chr}{$pos}{$ref}{$alt}{"format"}{$sample_id}{$format})){
							$sample_format .= ":".$h_result{$chr}{$pos}{$ref}{$alt}{"format"}{$sample_id}{$format};
						}
						else{
							$sample_format .= ":.";
						}
					}
					$sample_format =~ s/^://;
					print "\t".$sample_format;
				}
				print "\n";
			}
		}
	}
}
