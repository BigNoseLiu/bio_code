#!/usr/bin/env perl
use strict;
use Cwd qw(abs_path);
use File::Basename;
use Getopt::Long;

###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.0";
my $DATE = "2021-02-02";
###################################################################

my($DNA, $RNA, $index, $out, $version);
GetOptions(
	"DNA:s" => \$DNA,
	"RNA:s" => \$RNA,
	"index:s" => \$index,
	"out:s" => \$out,
	"version|v" => \$version,
);
if($version){
    print basename $0." $VERSION\n";
    exit(0);
}
&help unless($DNA && $RNA);
$index ||= "RPM";
$out ||= "./ActiveExpression.txt";
$out = abs_path $out;

my %threshold = (
	"Archaea" => 1.5,
	"Bacteria" => 1.5,
	"Eukaryota:Fungi" => 2,
	"Eukaryota:Parasite" => 2,
	"Eukaryota:Protozoa" => 2,
	"Viruses" => 1.5,
);

my %col = (
	"TPM" => 10,
	"RPM" => 11,
);

open IN,"$RNA" or die $!;
my %hash;
<IN>;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$hash{$t[5]} = $t[$col{$index}-1];
}
close IN;

open IN,"$DNA" or die $!;
open OU,">$out" or die $!;
<IN>;
print OU "Kingdom\tScientificName\tActiveness\n";
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if(exists $hash{$t[5]} && $hash{$t[5]}/$t[$col{$index}-1] >= $threshold{$t[0]}){
		print OU sprintf("$t[0]\t$t[5]\t%.2f\n", $hash{$t[5]}/$t[$col{$index}-1]);
	}
}
close IN;
close OU;

sub help{
print "
		Usage: perl $0

		--DNA        DNA sample's Detail.txt
		--RNA        RNA sample's Detail.txt
		--index      RPM or TPM [RPM]
		--out        [./ActiveExpression.txt]
		--version    Print version information.

";
exit(0);
}
