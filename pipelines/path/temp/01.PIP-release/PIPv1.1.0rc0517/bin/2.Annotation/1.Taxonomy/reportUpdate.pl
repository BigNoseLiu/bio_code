#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.1";
my $DATE = "2021-07-03";
###################################################################

my($raw, $update, $version);
GetOptions(
	"raw:s" => \$raw,
	"update:s" => \$update,
	"version|v" => \$version,
);

if($version){
    print basename $0." $VERSION\n";
    exit(0);
}

&help unless($raw && $update);
$raw = abs_path $raw;
$update = abs_path $update;

open IN,"$raw" or die $!;
my %hash;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$hash{$t[4]} = $t[3];
}
close IN;

open IN,"$update" or die $!;
open OU,">$update.tmp" or die $!;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if(exists $hash{$t[4]}){
		if($hash{$t[4]} eq 'S' && $t[3] eq '-'){
			$t[3] = 'S1';
		}else{
			$t[3] = $hash{$t[4]};
		}
		print OU (join "\t", @t)."\n";
	}else{
		print OU "$_\n";
	}
}
close IN;
`mv $update.tmp $update`;

sub help{
print "
	Usage: perl $0

	--raw                kraken's raw report file
	--update             kraken's filtered report file
	--version|v          print version information

";
exit(0);
}
