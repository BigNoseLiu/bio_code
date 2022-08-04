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
while( $line = <STDIN> ){
	chomp $line;

	#vcf headers
	if( $line =~ /^#/ ){
		$line =~ s/contig=<ID=chrM/contig=<ID=MT/;
		$line =~ s/contig=<ID=chr/contig=<ID=/;
		print $line."\n";
		next;
	}
	$line =~ s/^chrM/MT/i;
	$line =~ s/^chr//i;
	my @arr = split(/\t/,$line);
	#vcf format列的格式
	my @format_arr = split(/:/,$arr[8]);
	#获取GT所在序号
	my $gt_i = 10000;
	for( my$j=0;$j<scalar@format_arr;$j++ ){
		$gt_i = $j if($format_arr[$j] eq "GT");
	}
	#判断是否是multiple allele
	if($arr[4] =~ /,/){
		my %h_split = ();
		my @alt_arr = split(/,/,$arr[4]);
		#各allele的处理
		for( my$i=0;$i<scalar@alt_arr;$i++ ){
			#Alt列修改
			$h_split{$i}{"alt"} = $alt_arr[$i];
			foreach my $info ( split( /;/,$arr[7] ) ){
				if( $info =~ /^([^=]+)=(.*,.*)$/ ){
					my $title = $1;
					my $value = $2;
					my @value_arr = split(/,/,$value);
					$h_split{$i}{"info"} .= ";$title=$value_arr[$i]" if( defined( $value_arr[$i] ) );
				}
				else{
					$h_split{$i}{"info"} .= ";$info";
				}
			}

			my $x = $i+1;
			#各样检测结果处理
			for( my$m=9;$m<scalar@arr;$m++ ){
				#获取当前样本检测结果
				my @format_value_arr = split(/:/,$arr[$m]);

				#优先处理GT列，并初始化./.的其他列为.
				if( $gt_i != 10000 && $format_value_arr[$gt_i] =~ /^\d+/ ){
					#修改当前样本的GT的内容
					my %t_hash_geno  = ();
					my $t_count_geno = 0;
					foreach my $t_geno( split(/\//,$format_value_arr[$gt_i]) ){
						if(!defined( $t_hash_geno{$t_geno} )){
							$t_hash_geno{$t_geno}++;
							$t_count_geno++;
						}
					}
					if( defined($t_hash_geno{$x}) && $t_count_geno == 1 ){
						$format_value_arr[$gt_i] = "1/1";
					}
					elsif( defined($t_hash_geno{$x}) ){
						$format_value_arr[$gt_i] = "0/1";
					}
					else{
						$format_value_arr[$gt_i] = "./.";
					}
				}
				else{
					$format_value_arr[$gt_i] = "./.";
				}

				for( my$j=0;$j<scalar@format_arr;$j++ ){
					if( $format_value_arr[$gt_i] eq './.' ){
						$format_value_arr[$j] = ".";
						$format_value_arr[$gt_i] = './.' if( $j == $gt_i );
						$h_split{$i}{"format"}{$m} .= ":$format_value_arr[$j]";
					}
					elsif( $format_arr[$j] eq 'AD' ){
						my @t_arr = split(/,/,$format_value_arr[$j]);
						$h_split{$i}{"format"}{$m} .= ":$t_arr[0],$t_arr[$x]" if(defined($t_arr[$x]));
					}
					elsif( $format_arr[$j] eq 'AF' ){
						my @t_arr = split(/,/,$format_value_arr[$j]);
						$h_split{$i}{"format"}{$m} .= ":$t_arr[$i]" if(defined($t_arr[$i]));
					}
					else{
						$h_split{$i}{"format"}{$m} .= ":$format_value_arr[$j]";
					}
				}
			}
		}
		foreach my $i( sort {$a <=> $b} keys(%h_split) ){
			$arr[4] = $alt_arr[$i];
			$h_split{$i}{"info"} =~ s/^;//;
			$arr[7] = $h_split{$i}{"info"};
			for( my$m=9;$m<scalar@arr;$m++ ){
				$arr[$m] = $h_split{$i}{"format"}{$m};
				$arr[$m] =~ s/^://;
			}
			print join("\t",@arr)."\n";
		}
	}
	else{
		print $line."\n";
	}
}
