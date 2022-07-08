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
        Used to ...
Usage:  
        perl $0 [options]
        Options:
	     -help : reveal help info
	     -in  <str>  : input file path, format as below:
			column 1:	sample_id
	     -out  <str>  : output file path
        Example:
             perl $0 -out  output_file
Author & Contact:
	Mingming Liu
        mingming_liu\@shengtingroup.com
Last updated:
        2013-10-16
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	foreach my $t_cmd( @_ ){
		print STDERR $t_cmd."\n";
	}
	print STDERR color "reset";#change back the text color
}

my ( $help, $infile, $outfile, $clinvar_stat_file, $optional_cmd, $hpo_gene2disease_file);
GetOptions(
	"help"=>\$help,
	"in=s"=>\$infile,
	"out=s"=>\$outfile,
	"clinvar=s"=>\$clinvar_stat_file,
	"hpo=s"=>\$hpo_gene2disease_file,
);
open CLI,$clinvar_stat_file or die $!;
my $line;
my @lines = <CLI>;
my %h_clinvar_stat = ();
foreach $line(@lines){
	chomp $line;
	my @arr = split(/\t/,$line);
	next if( $arr[0] eq "-" );
	if( $arr[3] eq "mc1" ){
		$h_clinvar_stat{ $arr[1] }{ $arr[4] } = $arr[6].'|'.$arr[5];
	}
}
close CLI;

open HPO,$hpo_gene2disease_file or die $!;
@lines = <HPO>;
my %h_gene2disease = ();
foreach $line(@lines){
	chomp $line;
	next if($line =~ /^#/);
	my @arr = split(/\t/,$line);
	if( $arr[0] =~ /^OMIM/ ){
		$h_gene2disease{$arr[1]}{"OMIM"}{$arr[0]}++;
	}
	if( $arr[0] =~ /^ORPHA/ ){
		$h_gene2disease{$arr[1]}{"ORPHA"}{$arr[0]}++;
	}
}
close HPO;

@lines = <STDIN>;
foreach $line(@lines){

	if( $line =~ /^#/ ){
		if( $line =~ /^##INFO=<ID=ANN_STANDARD,/ ){
			print '##INFO=<ID=CLINVAR_STAT,Number=.,Type=String,Description="clinvar statistic info">'."\n";
			print '##INFO=<ID=OMIM,Number=.,Type=String,Description="OMIM diseases">'."\n";
			print '##INFO=<ID=ORPHA,Number=.,Type=String,Description="ORPHA diseases">'."\n";
		}
		print $line;
		next;
	}

	my $t_clinvar_stat = "-";
	my $t_omim = "-";
	my $t_orphan = "-";
	if( $line =~ /ANN_STANDARD=([^;\t]+)[;\s]/ ){
		my $t_snpeff_anno = $1;
		my @t_arr = split(/\|/,$t_snpeff_anno);
		if( defined( $h_clinvar_stat{$t_arr[6]}{$t_arr[1]} ) ){
			$t_clinvar_stat = $h_clinvar_stat{$t_arr[6]}{$t_arr[1]};
		}
		if( defined( $h_gene2disease{$t_arr[3]}{"OMIM"} ) ){
			$t_omim = join('|',keys(%{$h_gene2disease{$t_arr[3]}{"OMIM"}}));
		}
		if( defined( $h_gene2disease{$t_arr[3]}{"ORPHA"} ) ){
			$t_orphan = join('|',keys(%{$h_gene2disease{$t_arr[3]}{"ORPHA"}}));
		}
	}


	$line =~ s/ANN_STANDARD=/OMIM=$t_omim;ORPHA=$t_orphan;CLINVAR_STAT=$t_clinvar_stat;ANN_STANDARD=/;
	print $line;

}
