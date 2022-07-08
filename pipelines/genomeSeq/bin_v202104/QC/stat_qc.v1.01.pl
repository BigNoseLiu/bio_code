
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

	my ( $sample_id, $file_path ) = split(/:/,$data);
	if( defined($h_sample{$sample_id}) ){
		print STDERR "Err: dup sample $sample_id\n";
		exit(0);
	}
	$h_sample{$sample_id}++;
	open IN,$file_path or die $!;
	my $line;
	while( $line = <IN> ){
		next if( $line !~ /\S/ );
		chomp $line;
		if( $line =~ /number of reads = (\S.+)$/ ){
			$h_stat{"totalReads"}{$sample_id} = $1;
		}
		if( $line =~ /number of mapped reads = (\S.+)$/ ){
			$h_stat{"mappedReads"}{$sample_id} = $1;
		}
		if( $line =~ /number of duplicated reads \(flagged\) = (\S.+)$/ ){
			$h_stat{"duplReads"}{$sample_id} = $1;
		}
		if( $line =~ /median insert size = (\S.+)$/ ){
			$h_stat{"medianInsertSize"}{$sample_id} = $1;
		}
		if( $line =~ /mean coverageData = (\S.+)$/ ){
			$h_stat{"meanCoverage"}{$sample_id} = $1;
		}
		if( $line =~ /GC percentage = (\S.+)$/ ){
			$h_stat{"GCPerc"}{$sample_id} = $1;
		}
		if( $line =~ /There is a (\S+)% of reference with a coverageData >= 30X/ ){
			$h_stat{"10xCoverage"}{$sample_id} = $1;
		}
		if( $line =~ /There is a (\S+)% of reference with a coverageData >= 10X/ ){
			$h_stat{"30xCoverage"}{$sample_id} = $1;
		}
	}
	close IN;
}
my @headers = ("totalReads","mappedReads","duplReads","meanCoverage","10xCoverage","30xCoverage","medianInsertSize","GCPerc");
my @sample_ids  = sort {$a cmp $b} keys(%h_sample);
open OUT,">$out_file1" or die $!;
print OUT "#sample_id\t".join("\t",@headers)."\n";
foreach my $sample_id( @sample_ids ){
	print OUT $sample_id;
	foreach my $header( @headers ){
		if( defined($h_stat{$header}{$sample_id}) ){
			print OUT "\t".$h_stat{$header}{$sample_id};
		}
		else{
			print OUT "\t-";
		}
	}
	print  OUT "\n";
}
close OUT;
