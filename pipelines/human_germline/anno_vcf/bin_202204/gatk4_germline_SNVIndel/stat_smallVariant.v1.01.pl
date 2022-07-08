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

#input1: perl $0 sample1.filter.vcf sample2.filter.vcf .......  sampeN.filter.vcf >stat.xls
#input2: perl $0 id1:sample1.filter.vcf id2:sample2.filter.vcf .......  idN:sampeN.filter.vcf >stat.xls
my %h_sort_mut = ();
my %h_stat = ();
my %h_count= ();
my %h_tmb_stat = ();
my %h_sample = ();
my @types = ("other","het","hom");
my @info_titles = ("ANN_STANDARD","ExGeneSCombi","ExGeneSVar","gnomAD_genome_POPMAX","gnomAD_exome_POPMAX","civic-gene","civic-variant","civic-variant_civic_url:civic-rating:civic-evidence_level:civic-disease:civic-drugs:civic-evidence_type:civic-evidence_direction:civic-clinical_significance:civic-variant_origin","civicRegion4SmallMut_20200317","cosmic81","avsnp150","CLNALLELEID","CLNDN","CLNDISDB","CLNREVSTAT","CLNSIG","OMIM");
my $bed_file = shift@ARGV;
my $out_file1 = shift@ARGV;
my $out_file2 = shift@ARGV;
foreach my $data( @ARGV ){

	print STDERR $data."\n";
	my ( $sample_id, $call_type, $file_path ) = split(/:/,$data);
	next if(!-f $file_path);
	if( defined($h_sample{$sample_id}) ){
		print STDERR "Err: dup sample $sample_id\n";
		exit(0);
	}
	open IN,$file_path or die "no file $data\t$file_path\t$!";
	my $line;
	while( $line = <IN> ){
		chomp $line;
		next if( $line =~ /^#/ );
		my @arr = split(/\t/,$line);
		next if( !defined($arr[9]) );
		next if( !defined($arr[6]) );
		#next if( $arr[6] !~ /PASS/i );
		#统一染色体命名法
		$line =~ s/^chr//i;
		my $type = $types[0];
		if( $arr[9] =~ /^0\/1:/i || $arr[9] =~ /^1\/0:/i ){
			$type = $types[1];
		}
		elsif( $arr[9] =~ /^1\/1:/i ){
			$type = $types[2];
		}
		if( $call_type =~ /somatic/i ){
			my %h_format_ids = ();
			my $format_id_count = 0;
			foreach my $format_id( split(/:/,$arr[8]) ){
				$h_format_ids{$format_id} = $format_id_count;
				$format_id_count++;
			}
			my @datas =  split(/:/,$arr[9]);
			my $dp = $datas[$h_format_ids{"DP"}];
			my @ads = split(/,/,$datas[$h_format_ids{"AD"}]);
			my $ad = $ads[-1];
			my $af = ($ad/$dp)*100;
			$af =~ s/(\.\d\d)\d*$/$1/;
			$type = $af."/".$dp;
		}
		if( !defined($h_count{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{"info_detail"}) ){
			$h_count{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{"info_detail"} = $arr[7];
			foreach my $info_title( @info_titles ){
				$h_count{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{"info"} .= "\t";
				if( $arr[7] =~ /$info_title=([^;]+);/ ){
					$h_count{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{"info"} .= "$1";
				}
			}
		}
		#过滤
		if( $arr[6] !~ /PASS/i ){
			$type = "Filtered";
			$h_count{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{"Filtered"}{$sample_id}++;
		}
		else{
			$h_count{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{"all"}{$sample_id}++;
			#统计somatic突变个数
			$h_sample{$sample_id}{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}++;
			#if( $call_type =~ /somatic/i ){
			#	$h_sample{$sample_id}{"somatic"}{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}++;
			#}
			#else{
			#	$h_sample{$sample_id}{"germline"}{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}++;
			#}
		}
		my $num_chr = $arr[0];
		$num_chr = 23 if($arr[0] =~ /x/i);
		$num_chr = 24 if($arr[0] =~ /y/i);
		$h_sort_mut{$num_chr}{$arr[1]}{$arr[3]}{$arr[4]} = join("\t",$arr[0],$arr[1],$arr[3],$arr[4]);
		$h_stat{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{$sample_id}{"genotype"} = $type;
		$h_stat{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{$sample_id}{"filter"} = $arr[6];
		$h_stat{join("\t",$arr[0],$arr[1],$arr[3],$arr[4])}{$sample_id}{"info"} = $arr[7];
	}
	close IN;
}
#my @infos = ( "genotype","filter","info" );
my @infos = ( "genotype","filter" );
my @sample_ids = sort {$a cmp $b} keys(%h_sample);
open OUT,">$out_file1" or die $!;
print OUT "#chr\tpos\tref\tmut\tcount_all\tall_sample\tcount_filter\tfilter_sample\t".join("\t",@info_titles);
foreach my $info( @infos ){
	foreach my $sample_id( @sample_ids ){
		if( $info eq "genotype" ){
			print OUT "\t$sample_id";
		}
		else{
			print OUT "\t$info\_$sample_id";
		}
	}
}
print OUT "\tdetail_info\n";
foreach my $chr( sort {$a <=> $b} keys(%h_sort_mut) ){
	foreach my $pos( sort {$a <=> $b} keys(%{$h_sort_mut{$chr}}) ){
		foreach my $ref( sort {$a cmp $b} keys(%{$h_sort_mut{$chr}{$pos}}) ){
			foreach my $alt( sort {$a cmp $b} keys(%{$h_sort_mut{$chr}{$pos}{$ref}}) ){

				my $mut = $h_sort_mut{$chr}{$pos}{$ref}{$alt};#join("\t",$chr,$pos,$ref,$alt);
				my @count_all_samples =  sort {$a cmp $b} keys(%{$h_count{$mut}{"all"}});
				my @count_filtered_samples = ();
				if( defined($h_count{$mut}{"Filtered"}) ){
					@count_filtered_samples = sort {$a cmp $b} keys(%{$h_count{$mut}{"Filtered"}});
				}
				print OUT $mut."\t".(scalar@count_all_samples)."\t".join(",",@count_all_samples)."\t".(scalar@count_filtered_samples)."\t".join(",",@count_filtered_samples).$h_count{$mut}{"info"};
				foreach my $info( @infos ){
					foreach my $sample_id( @sample_ids ){
						if( defined($h_stat{$mut}{$sample_id}{$info}) ){
							print OUT "\t".$h_stat{$mut}{$sample_id}{$info};
						}
						else{
							print OUT "\t-";
						}
					}
				}
				print OUT "\t".$h_count{$mut}{"info_detail"}."\n";
			}


		}
	}
}
close OUT;

my %h_bed = ();
open BED, $bed_file or die $!;
my @bed_lines = <BED>;
foreach my $line( @bed_lines ){
	chomp $line;
	$line =~ s/^chr//i;
	my @arr = split(/\t/,$line);
	for( my $i=$arr[1]+1;$i<=$arr[2];$i++ ){
		$h_bed{$arr[0]."\t".$i}++;
	}
}
close BED;
my $bed_count = scalar(keys(%h_bed));
print "bed count:\t$bed_count\n";
open OUT,">".$out_file2 or die $!;
print OUT "#sample\ttarget_region\tmut_count\tTMB\n";
foreach my $sample_id( @sample_ids ){
	#if( defined($h_sample{$sample_id}{"somatic"}) ){
		#my @arr_variant = sort {$a cmp $b} keys(%{$h_sample{$sample_id}{"somatic"}});
		my @arr_variant = sort {$a cmp $b} keys(%{$h_sample{$sample_id}});
		my $mut_count = scalar@arr_variant;
		my $TMB = $mut_count/$bed_count*1000000;
		$TMB =~ s/(\.\d\d)\d*$/$1/;
		print OUT $sample_id."\t$bed_count\t$mut_count\t$TMB\n";
	#}
}
close OUT;
