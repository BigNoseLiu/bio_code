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
        Used to merge all annotation infomation.
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

my ($help,$exomiser_file,$ref_fa_file,$outfile);
GetOptions(
	"exomiser=s"=>\$exomiser_file,
	"help"=>\$help,
);
if (defined $help ) {
	&usageWithColor();
	exit 0;
}

my $line;
#加载exomiser结果
my %h_exomser = ();
my @exomiser_header = ();
if( defined($exomiser_file) ){
	open IN, $exomiser_file or die $!;
	$line = <IN>;	chomp $line;
	my @arr = split(/\t/,$line);
	print join("\t",$arr[28],$arr[29],$arr[30],$arr[31])."\t";
	while( $line = <IN> ){
		chomp $line;
		@arr = split(/\t/,$line);
		$h_exomser{join("\t",$arr[0],$arr[1],$arr[2],$arr[3])} = join("\t",$arr[28],$arr[29],$arr[30],$arr[31]);
	}
	close IN;
}

#加载annovar.vcf结果，获取纯合、杂合类型
my %h_vcf = ();
my %h_sample = ();
my $annovar_vcf = $ARGV[0] or die $!;
my @vcf_header = ();
open IN, $annovar_vcf or die $!;
while( $line = <IN> ){
	chomp $line;
	#vcf headers
	if( $line =~ /^#/ ){
		if( $line =~ /^#CHROM/ ){
			@vcf_header = split( /\t/, $line );
		}
		next;
	}
	my @arr = split(/\t/,$line);
	#vcf format列的格式
	my @format_arr = split(/:/,$arr[8]);
	#获取GT所在序号
	my $gt_i = 10000;
	my $af_i = 10000;
	for( my$j=0;$j<scalar@format_arr;$j++ ){
		$gt_i = $j if($format_arr[$j] eq "GT");
		$af_i = $j if($format_arr[$j] eq "AF");
	}
	for( my$m=9;$m<scalar@arr;$m++ ){
		#获取当前样本检测结果
		my @format_value_arr = split(/:/,$arr[$m]);

		my $mut_type = "wild";
		#优先处理GT列，并初始化./.的其他列为.
		if( $gt_i != 10000 && $format_value_arr[$gt_i] =~ /^0\/1/ ){
			$mut_type = "het";
		}
		elsif( $gt_i != 10000 && $format_value_arr[$gt_i] =~ /^1\/1/ ){
			$mut_type = "hom";
		}
		if( $arr[0] =~ /^chrm/i || $arr[0] =~ /^m/i ){
			if( $af_i != 10000 && $format_value_arr[$af_i] =~ /\d/ ){
				$mut_type = $format_value_arr[$af_i];
			}
		}
		$h_vcf{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{$vcf_header[$m]} = join("\t",$mut_type,$arr[$m]);
		$h_sample{$vcf_header[$m]}++;
	}
}
close IN;
foreach my $sample( sort {$a cmp $b} keys(%h_sample) ){
	print "mut_result_$sample\tmut_detail_$sample\t";
}

#加载annovar.txt结果，获取纯合、杂合类型
my $annovar_txt = $ARGV[1] or die $!;
open IN, $annovar_txt or die $!;

#获取标题所对应的列序号
$line = <IN>;	chomp $line;
my @arr_header = split(/\t/,$line);

#修改校正annovar的标题
my %h_count = ();
for(my$i=0;$i<@arr_header;$i++){
	my $t_header = $arr_header[$i];
	$t_header =~ s/_SCORE/_score/;
	$t_header = "GERP++_score" if($t_header eq "GERP++_RS");
	$t_header = "CADD_score" if($t_header eq "CADD_raw");
	$t_header = "Eigen_score" if($t_header eq "Eigen-PC-raw_coding");
	$t_header =~ s/^AF$/gnomadGenome_all/;
	$t_header =~ s/^AF_(\S+)$/gnomadGenome_$1/;
	if( defined($h_count{$t_header}) ){
		$h_count{$t_header}++;
		$t_header .= "-".$h_count{$t_header};
		$t_header =~ s/^gnomadGenome_(\S+)-2$/gnomadExome_$1/;
	}
	$arr_header[$i] = $t_header;
	$h_count{$t_header}++;
}

print join("\t",@arr_header)."\n";

my %h_header = ();
for(my$i=0;$i<@arr_header;$i++){
	my $t_header = $arr_header[$i];
	$h_header{$t_header} = $i;
}



while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	my $mut_id = join("\t",$arr[$h_header{'Otherinfo4'}],$arr[$h_header{'Otherinfo5'}],$arr[$h_header{'Otherinfo7'}],$arr[$h_header{'Otherinfo8'}]);
	#注释exomiser结果
	if( defined($exomiser_file) ){
		my $exomiser_result = join("\t","-1","-1","-1","-1");
		if( defined($h_exomser{$mut_id}) ){
			$exomiser_result = $h_exomser{$mut_id};
		}
		print $exomiser_result."\t";
	}
	#注释各样本检出结果
	foreach my $sample( sort {$a cmp $b} keys(%h_sample) ){
		my $mut_result = join("\t","-","-");
		if( defined($h_vcf{$mut_id}) ){
			$mut_result = $h_vcf{$mut_id}{$sample};
		}
		print "$mut_result\t";
	}
	print $line."\n";
}

close IN;
