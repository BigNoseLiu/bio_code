#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);



my $line;
while( $line = <STDIN> ){
	next if($line =~ /^##/);
	chomp $line;
	my @arr = split(/\t/,$line);
	$arr[7] =~ s/^INFO$/ExonicFunc\tAAChange/;
	$arr[7] =~ s/^.*ExonicFunc.refGeneWithVer=([^;]*);AAChange.refGeneWithVer=([^;]*);.*$/$1\t$2/;
	print join("\t",@arr)."\n";
}
