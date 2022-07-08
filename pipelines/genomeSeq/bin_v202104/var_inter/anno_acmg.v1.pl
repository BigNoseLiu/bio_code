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

my ( $help, $infile, $out_dir, $GATK_PON, $optional_cmd, $debug );
GetOptions(
	"help"=>\$help,
	"debug"=>\$debug,
	"opt"=>\$optional_cmd,
	"in=s"=>\$infile,
	"out=s"=>\$out_dir,
	"pon=s"=>\$GATK_PON,
);
my @acmg_evidence = ("PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7");
my %h_var_acmg = ();
