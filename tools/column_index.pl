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
        Used to get the column number for each field of one line.
Usage:  
        perl $0 [options]
        Options:
	     -help : reveal help info
	     -sep  <str>  : field seperator. Default: \\t
        Example:
             perl $0 <in_file -sep " "
Author & Contact:
	Mingming Liu
Last updated:
        2016-12-6
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}
my %h_transfer= (
	'\t'=>"\t",
	'\\\\'=>'\\',
	'\n'=>"\n"
);

my ($help);
my $seperator = "\t";
GetOptions(
	"help"=>\$help,
	"sep=s"=>\$seperator,
);
$seperator =~ s/(\\.)/$h_transfer{$1}/g;


if (defined $help){
	&usageWithColor();
	exit 0;
}

#read header
my $line = <STDIN>;
chomp $line;
my @arr1 = ();
@arr1 = split($seperator,$line) if(defined( $line ));

#read content
my @arr2 = ();
$line = <STDIN>;
chomp $line;
if(defined( $line )){
	@arr2 = split($seperator,$line);
	chomp $line;
}
my $count = 1;
for( my $i=0;$i<scalar @arr1 || $i<scalar @arr2;$i++ ){
	$arr1[$i] = " " if( !defined( $arr1[$i] ) );
	print $i."\t".$arr1[$i];
	$arr2[$i] = " " if( !defined( $arr2[$i] ) );
	print "\t:\t".$arr2[$i]."\n";
}
