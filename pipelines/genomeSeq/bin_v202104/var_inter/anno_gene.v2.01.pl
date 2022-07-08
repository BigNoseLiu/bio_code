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
        2020-03-30
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

my ( $help, $infile, $outfile, $clinvar_stat_file, $clinvar_stat_multiple_file, $optional_cmd, $hpo_gene2phenotype_file, $hpo_phenotype2gene_file, $exomiser_file);
GetOptions(
	"help"=>\$help,
	"in=s"=>\$infile,
	"out=s"=>\$outfile,
	"clinvar_all=s"=>\$clinvar_stat_file,
	"exomiser=s"=>\$exomiser_file,
	"clinvar_mul=s"=>\$clinvar_stat_multiple_file,
	"hpo_g2p=s"=>\$hpo_gene2phenotype_file,
	"hpo_p2g=s"=>\$hpo_phenotype2gene_file,
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
		$h_clinvar_stat{ $arr[1] }{"all"}{ $arr[4] } = $arr[6].'|'.$arr[5];
	}
}
close CLI;


open CLI,$clinvar_stat_multiple_file or die $!;
@lines = <CLI>;
foreach $line(@lines){
	chomp $line;
	my @arr = split(/\t/,$line);
	next if( $arr[0] eq "-" );
	if( $arr[3] eq "mc1" ){
		$h_clinvar_stat{ $arr[1] }{"mul"}{ $arr[4] } = $arr[6].'|'.$arr[5];
	}
}
close CLI;


my %h_exomiser = ();
my $exomiser_header_vcf = "";
if(defined($exomiser_file) && $exomiser_file =~ /\S/){
	open CLI,$exomiser_file or die $!;
	@lines = <CLI>;
	my @to_anno_exomser_headers = ("ExGeneSCombi","ExGeneSPheno","ExGeneSVar","ExVarScore");
	foreach $line(@lines){
		chomp $line;
		if( $line =~ /^#/ ){
			foreach my $header( @to_anno_exomser_headers ){
				if( $line =~ /$header/ ){
					$exomiser_header_vcf .= "$line\n";
				}
			}
			next;
		}
		$line =~ s/^chr//i;
		my @arr = split(/\t/,$line);
		my $mut = join("\t",$arr[0],$arr[1],$arr[3],$arr[4]);
		$h_exomiser{ $mut } = "";
		foreach my $header( @to_anno_exomser_headers ){
			if( $line =~ /($header=[^;]+);/ ){
				$h_exomiser{ $mut } .= "$1;";
			}
		}
	}
	close CLI;
}

my %h_gene2hpo = ();
my %h_gene2disease = ();
open HPO,$hpo_phenotype2gene_file or die $!;
@lines = <HPO>;
foreach $line(@lines){
	chomp $line;
	next if($line =~ /^#/);
	my @arr = split(/\t/,$line);
	if( $arr[6] =~ /^OMIM/ ){
		$h_gene2disease{$arr[3]}{"OMIM"}{$arr[6]}++;
	}
	if( $arr[6] =~ /^ORPHA/ ){
		$h_gene2disease{$arr[3]}{"ORPHA"}{$arr[6]}++;
	}
	$h_gene2hpo{$arr[3]}{"all"}{$arr[0]}++;
}
close HPO;


open HPO,$hpo_gene2phenotype_file or die $!;
@lines = <HPO>;
foreach $line(@lines){
	chomp $line;
	next if($line =~ /^#/);
	my @arr = split(/\t/,$line);
	if( $arr[8] =~ /^OMIM/ ){
		$h_gene2disease{$arr[1]}{"OMIM"}{$arr[8]}++;
	}
	if( $arr[8] =~ /^ORPHA/ ){
		$h_gene2disease{$arr[1]}{"ORPHA"}{$arr[8]}++;
	}
	$h_gene2hpo{$arr[1]}{"all"}{$arr[2]}++;
	$h_gene2hpo{$arr[1]}{"dir"}{$arr[2]}++;
}
close HPO;


@lines = <STDIN>;
foreach $line(@lines){

	if( $line =~ /^#/ ){
		if( $line =~ /^##INFO=<ID=ANN_STANDARD,/ ){
			print '##INFO=<ID=CLINVAR_STAT_ALL,Number=.,Type=String,Description="clinvar statistic info for all variants">'."\n";
			print '##INFO=<ID=CLINVAR_STAT_MUL,Number=.,Type=String,Description="clinvar statistic info for variants with multiple support">'."\n";
			print '##INFO=<ID=OMIM,Number=.,Type=String,Description="OMIM diseases">'."\n";
			print '##INFO=<ID=ORPHA,Number=.,Type=String,Description="ORPHA diseases">'."\n";
			print '##INFO=<ID=HPO_DIR,Number=.,Type=String,Description="HPO_phenotype directly related with this gene">'."\n";
			print '##INFO=<ID=HPO_ALL,Number=.,Type=String,Description="HPO_phenotype related with this gene(including the parent hpo item)">'."\n";
			if( $exomiser_header_vcf =~ /\S/ ){
				print $exomiser_header_vcf;
			}
		}
		print $line;
		next;
	}
	$line =~ s/^chr//i;
	my @arr = split(/\t/,$line);
	my $mut = join("\t",$arr[0],$arr[1],$arr[3],$arr[4]);

	my $t_clinvar_stat_all = "-";
	my $t_clinvar_stat_mul = "-";
	my $t_omim = "-";
	my $t_orphan = "-";
	my $t_hpo_all = "-";
	my $t_hpo_dir = "-";
	if( $line =~ /ANN_STANDARD=([^;\t]+)[;\s]/ ){
		my $t_snpeff_anno = $1;
		my @t_arr = split(/\|/,$t_snpeff_anno);
		if( defined( $h_clinvar_stat{$t_arr[6]}{"all"}{$t_arr[1]} ) ){
			$t_clinvar_stat_all = $h_clinvar_stat{$t_arr[6]}{"all"}{$t_arr[1]};
		}
		if( defined( $h_clinvar_stat{$t_arr[6]}{"mul"}{$t_arr[1]} ) ){
			$t_clinvar_stat_mul = $h_clinvar_stat{$t_arr[6]}{"mul"}{$t_arr[1]};
		}
		if( defined( $h_gene2disease{$t_arr[3]}{"OMIM"} ) ){
			$t_omim = join('|',keys(%{$h_gene2disease{$t_arr[3]}{"OMIM"}}));
		}
		if( defined( $h_gene2disease{$t_arr[3]}{"ORPHA"} ) ){
			$t_orphan = join('|',keys(%{$h_gene2disease{$t_arr[3]}{"ORPHA"}}));
		}
		if( defined( $h_gene2hpo{$t_arr[3]}{"all"} ) ){
			$t_hpo_all = join('|',keys(%{$h_gene2hpo{$t_arr[3]}{"all"}}));
		}
		if( defined( $h_gene2hpo{$t_arr[3]}{"dir"} ) ){
			$t_hpo_dir = join('|',keys(%{$h_gene2hpo{$t_arr[3]}{"dir"}}));
		}
	}


	if( defined( $h_exomiser{ $mut } ) && $h_exomiser{$mut} =~ /\S/ ){
		my $exomiser_info = $h_exomiser{ $mut };
		$exomiser_info =~ s/;$//;
		$line =~ s/ANN_STANDARD=/$exomiser_info;ANN_STANDARD=/;
	}
	$line =~ s/ANN_STANDARD=/OMIM=$t_omim;ORPHA=$t_orphan;CLINVAR_STAT_ALL=$t_clinvar_stat_all;CLINVAR_STAT_MUL=$t_clinvar_stat_mul;ANN_STANDARD=/;
	$line =~ s/(ANN_STANDARD=[^\t]+);?\t/$1;HPO_DIR=$t_hpo_dir;HPO_ALL=$t_hpo_all\t/;
	print $line;

}
