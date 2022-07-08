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
my $line = $ARGV[0];
while( $line =~ s/(.)$// ){
	print $h_reverse{$1};
}
print "\n";
