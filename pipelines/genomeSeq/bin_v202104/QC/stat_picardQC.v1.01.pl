
#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %h_stat = ();
my %h_sample = ();
my $out_file1 = shift@ARGV;
foreach my $data( @ARGV ){

	my @input_files = split(/:/,$data);
	my $sample_id = shift @input_files; 
	my $file_path = join(" ",@input_files);

	$h_sample{$sample_id}++;
	open IN,"cat $file_path|" or die $!;
	my $line;
	while( $line = <IN> ){
		chomp $line;



		if( $line =~ /^##\s*METRICS CLASS\s*picard.analysis.AlignmentSummaryMetrics/i ){
			my $t_line = <IN>;
			my @headers = split(/\t/,$t_line);
			$t_line = <IN>;
			$t_line = <IN>;	chomp $t_line;
			my @t_arr = split(/\t/,$t_line);
			for( my $i=0;$i<scalar@headers;$i++ ){
				$h_stat{$headers[$i]}{$sample_id} = $t_arr[$i];
			}
		}

		if( $line =~ /^##\s*METRICS CLASS\s*picard.sam.DuplicationMetrics/i ){
			my $t_line = <IN>;
			my @headers = split(/\t/,$t_line);
			$t_line = <IN>;	chomp $t_line;
			my @t_arr = split(/\t/,$t_line);
			for( my $i=0;$i<scalar@headers;$i++ ){
				$h_stat{$headers[$i]}{$sample_id} = $t_arr[$i];
			}
		}



		if( $line =~ /^##\s*METRICS CLASS\s*picard.analysis.InsertSizeMetrics/i ){
			my $t_line = <IN>;
			my @headers = split(/\t/,$t_line);
			$t_line = <IN>;	chomp $t_line;
			my @t_arr = split(/\t/,$t_line);
			for( my $i=0;$i<scalar@headers;$i++ ){
				$h_stat{$headers[$i]}{$sample_id} = $t_arr[$i];
			}
		}

		if( $line =~ /^##\s*METRICS CLASS\s*picard.analysis.CollectQualityYieldMetrics\$QualityYieldMetrics/i ){
			my $t_line = <IN>;
			my @headers = split(/\t/,$t_line);
			$t_line = <IN>;	chomp $t_line;
			my @t_arr = split(/\t/,$t_line);
			for( my $i=0;$i<scalar@headers;$i++ ){
				$h_stat{$headers[$i]}{$sample_id} = $t_arr[$i];
			}
		}
		if( $line =~ /^##\s*METRICS CLASS\s*picard.analysis.directed.HsMetrics/i ){
			my $t_line = <IN>;
			my @headers = split(/\t/,$t_line);
			$t_line = <IN>;	chomp $t_line;
			my @t_arr = split(/\t/,$t_line);
			for( my $i=0;$i<scalar@headers;$i++ ){
				$h_stat{$headers[$i]}{$sample_id} = $t_arr[$i];
			}
		}


	}
	close IN;
}
#my @headers = ("totalReads","mappedReads","duplReads","meanCoverage","10xCoverage","30xCoverage","medianInsertSize","GCPerc");
my @headers = ("TOTAL_READS","READ_LENGTH","TOTAL_BASES","PCT_PF_READS_ALIGNED","Q30_BASES","MEDIAN_INSERT_SIZE","PERCENT_DUPLICATION","TARGET_TERRITORY","PCT_SELECTED_BASES","MEAN_TARGET_COVERAGE","PCT_TARGET_BASES_1X","PCT_TARGET_BASES_2X","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_100X");
my @sample_ids  = sort {$a cmp $b} keys(%h_sample);
open OUT,">$out_file1" or die $!;
my $out_header = join("\t",@headers);
$out_header =~ s/SELECTED_BASES/CaptureRatio/g;
$out_header =~ s/PCT_PF_/\(%\)/g;
$out_header =~ s/PCT_/\(%\)/g;
$out_header =~ s/PERCENT_/\(%\)/g;
$out_header =~ s/Q30_BASES/\(%\)Q30/g;
print OUT "#sample_id\t$out_header\n";
foreach my $sample_id( @sample_ids ){
	print OUT $sample_id;
	foreach my $header( @headers ){
		if( defined($h_stat{$header}{$sample_id}) ){
			if( $header eq "MEAN_TARGET_COVERAGE"){
				$h_stat{$header}{$sample_id} =~ s/(\.\d\d).*$/$1/;
			}
			if( $header =~ /PCT_TARGET_BASES/ || $header eq "PCT_PF_READS_ALIGNED" || $header eq "PERCENT_DUPLICATION"){
				$h_stat{$header}{$sample_id} *= 100;
				$h_stat{$header}{$sample_id} =~ s/(\.\d\d).*$/$1/;
			}
			if( $header eq "Q30_BASES" ){
				$h_stat{$header}{$sample_id} = $h_stat{$header}{$sample_id}/$h_stat{"TOTAL_BASES"}{$sample_id}*100;
				$h_stat{$header}{$sample_id} =~ s/(\.\d\d).*$/$1/;
			}
			print OUT "\t".$h_stat{$header}{$sample_id};
		}
		else{
			print OUT "\t-";
		}
	}
	print  OUT "\n";
}
close OUT;
