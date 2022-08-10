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
        Used to locate the positon of the seq in reference sequence.
Reference genomes:
	/ifs1/NGB_DB/liumingming/XArchive4Gene/database/ftp_oriented/virus/All/virus.all.fa
	/ifs1/NGB_DB/liumingming/XArchive4Gene/database/ftp_oriented/genome/human/hg19/hg19.fa
Usage:  
        perl $0 [options]
        Options:
	     -help : reveal help info
        Example:
             echo ATCGG GGATC |perl $0 -ref hg19.fa
Author & Contact:
	Mingming Liu
        liumingming\@genomics.cn
Last updated:
        2014-2-27
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}

my ($help,$not_fuzzy_match,$number_chr,$ref_files,$chr,$pos);
#default settings
$ref_files = "/ifs1/NGB_DB/liumingming/XArchive4Gene/database/ftp_oriented/genome/human/hg19/hg19.fa";
GetOptions(
	"help"=>\$help,
	"whole"=>\$not_fuzzy_match,
	"number"=>\$number_chr,
	"ref=s"=>\$ref_files,
	"chr=s"=>\$chr,
	"pos=s"=>\$pos,
);
if (defined $help){
	&usageWithColor();
	exit 0;
}
my $line;

my %h_chr_fa = ();
foreach my $t_ref( split(/,/,$ref_files) ){
	open REF,"$t_ref" or die "Fail opening $t_ref\n";
	$/ = '>';
	my $title = "";
	my $chr= 0;
	while( $line = <REF> ){
		chomp $line;
		next if( $line =~ /^\s*$/ );
		if( $line =~ s/^([^\n]*)\n// ){
			$title = $1;
			$chr = $title;
			if( $title =~ /^(\S+)\s/ ){
				$chr = $1;
			}
		}
		else{
			print STDERR "wrong formated fa:$t_ref\n";
			next;
		}
		$line =~ s/[>\s]//g;
		$h_chr_fa{$chr} = $line;
	}
	close REF;
}

my %h_reverse = (
	"A"=>"T",
	"a"=>"t",
	"T"=>"A",
	"t"=>"a",
	"C"=>"G",
	"c"=>"g",
	"G"=>"C",
	"g"=>"c",
);



while( $line = <STDIN> ){
	chomp $line;
	foreach my $seq( split(/\s+/,$line) ){
		next if( $seq !~ /\S/ );
		my $locate = "-";
		my $locate_count = 0;

		my $temp_seq = $seq;
		my $rev_seq = "";
		my $rev_locate = "-";
		my $rev_locate_count = 0;
		while( $temp_seq =~ s/(.)$// ){
			my $t_base = $1;
			if( !defined($h_reverse{$t_base}) ){
				$h_reverse{$t_base} = $t_base;
			}
			$rev_seq .= $h_reverse{$t_base};
		}

		foreach my $chr( sort {$a cmp $b} keys(%h_chr_fa)){
			while( $h_chr_fa{$chr} =~ /($seq)/ig ){
				$locate .= ",$chr:".($-[0]+1)."_".($-[0]+length($seq)) if( $locate_count < 3 );
				$locate_count++;
			}
		}


		foreach my $chr( sort {$a cmp $b} keys(%h_chr_fa)){
			while( $h_chr_fa{$chr} =~ /($rev_seq)/ig ){
				$rev_locate .= ",$chr:".($-[0]+1)."_".($-[0]+length($rev_seq)) if( $rev_locate_count < 3 );
				$rev_locate_count++;
			}
		}

		$locate =~ s/^-,//;
		$rev_locate =~ s/^-,//;
		print $seq."\t$locate\t$rev_locate\n";
	}
}
