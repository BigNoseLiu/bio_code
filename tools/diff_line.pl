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
        Used to check if the field values of different files match with the first file;
Usage:  
        perl $0 column1:file1 column2:file2 ....
        Options:
	     -help : reveal help info
	     column:file
	     	column: the 1_based column number of the field you want to check
	     	file: file name
        Example:
             perl $0 1:file1 3:file2 7:file3
Author & Contact:
	Mingming Liu
        liumingming\@genomics.cn
Last updated:
        2013-10-16
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}

my ($help,$out);
GetOptions(
	"help"=>\$help,
	"out=s"=>\$out,
);
if (defined $help || scalar @ARGV < 2 ) {
	&usageWithColor();
	exit 0;
}

#print "Start time:".time."\n";
my $file_count = 0;
my %h_first;
my $first_file_name = "";
foreach my $t_record( @ARGV ){
	my @a_record= split( /:/,$t_record );
	if( scalar( @a_record ) != 2 || ($a_record[0] =~ /\D/ || $a_record[0] =~ /^$/) || !(-f $a_record[1])){
		&usageWithColor();
		exit 0;
	}
	$file_count++;
	my ($i,$file_name) = @a_record;
	$i--;
	open IN,"$file_name" or die "Fail opening $file_name\n";
	my $line;
	while( $line = <IN> ){
		chomp $line;
		my @arr = split(/\t/,$line);
		$arr[$i] =~ s/^\s*(\S)/$1/;
		$arr[$i] =~ s/(\S)\s*$/$1/;
		$h_first{$arr[$i]}{$t_record}++;
	}
	close IN;
}
foreach my $t_fieldValue( keys( %h_first ) ){
	my @t_files = keys( %{$h_first{$t_fieldValue}} );
	if( scalar @t_files == 1 ){
		print $t_fieldValue."\tis unique in\t".$t_files[0]."\n";
	}
}
