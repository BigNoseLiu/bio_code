#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);

sub usageWithColor {
	print color "green";#change the text color
	print <<USAGE;
Description:
        Used to find out the alias names of human gene.
	Gene Informatioin Source:
		ftp.ncbi.nih.gov/gene/DATA/20131211/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
Usage:  perl $0 [-nc]  in_names
        Options:
	     -help : reveal help info
	     -nc : if defined, No color will be added
	     in_names: names of human gene.(seperated by ':')
        Example:
             perl $0 AQP1:CO:AADAC
             perl $0 -nc AQP1:CO:AADAC:Liumm
Author & Contact:
	Mingming Liu
        liumingming\@genomics.cn
Last updated:
        2014-2-12
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}

my ($help,$no_color);
GetOptions(
	"help"=>\$help,
	"nc"=>\$no_color,
);
if (defined $help || scalar @ARGV < 1 ) {
	&usageWithColor();
	exit 0;
}

my %h_geneList;
foreach my $geneName( split(/:/,$ARGV[0]) ){
	$h_geneList{ $geneName } = 0;
}

my $gene_info_file = "$Bin/../database/ftp_oriented/gene_info_NCBI/Homo_sapiens.gene_info";
my $line;
my %h_geneInfo;
open GENE,"$gene_info_file" or die "Fail opening $gene_info_file\n";
while( $line = <GENE> ){
	chomp $line;
	next if( $line =~ /^#/ );
	my @arr = split( /\t/,$line );
	my @arr4 = split(/\|/,$arr[4]);
	my $flag = 0;
	if( defined( $h_geneList{ $arr[2] } ) ){
		$flag = 1;
	}
	else{
		foreach my $gene_name( @arr4 ){
			if( defined( $h_geneList{ $gene_name } ) ){
				$flag = 1;
				last;
			}
		}
	}
	if( $flag == 1 ){
		foreach my $gene_name( @arr4 ){
			$h_geneInfo{ $arr[2] }{ $gene_name }++;
		}
	}
}
close GENE;


sub print_red{
	if( !defined( $no_color ) ){
		print color "red";#change the text color
		print $_[0];
		print color "reset";#change back the text color
	}
	else{
		print $_[0];
	}
}


sub print_green{
	if( !defined( $no_color ) ){
		print color "green";#change the text color
		print $_[0];
		print color "reset";#change back the text color
	}
	else{
		print $_[0];
	}
}

foreach my $gene_name( sort {$a cmp $b} keys( %h_geneInfo ) ){
	if( defined( $h_geneList{ $gene_name } ) ){
		&print_green( $gene_name );
		$h_geneList{ $gene_name } = 1;
	}
	else{
		print $gene_name;
	}
	print "\t";
	my $flag = 0;
	foreach my $gene_name2( sort {$a cmp $b} keys( %{ $h_geneInfo{ $gene_name } } ) ){
		if( $flag == 1 ){
			print "|";
		}
		$flag = 1;
		if( defined( $h_geneList{ $gene_name2 } ) ){
			&print_green( $gene_name2 );
			$h_geneList{ $gene_name2 } = 1;
		}
		else{
			print $gene_name2;
		}
	}
	print "\n";
}

foreach my $gene_name( sort {$a cmp $b} keys( %h_geneList) ){
	if( $h_geneList{ $gene_name } == 0 ){
		&print_red( $gene_name."\t-NA-\n" );
	}
}
