
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
	my $line = <IN>;
	chomp $line;
	my @headers = split(/\t/,$line);
	while( $line = <IN> ){
		chomp $line;
		next if( $line !~ /\S/ );
		my @arr = split(/\t/,$line);
		for( my$i=0;$i<scalar@headers;$i++ ){
			$h_stat{$headers[$i]}{$sample_id} = $arr[$i];
		}
	}
	close IN;
}
my @headers = sort {$a cmp $b} keys(%h_stat);
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
