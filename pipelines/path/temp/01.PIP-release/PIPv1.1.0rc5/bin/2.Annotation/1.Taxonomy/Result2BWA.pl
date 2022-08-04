#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../../../lib/perl5";
use MyModule::GlobalVar qw($PIGZ);
###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.0";
my $DATE = "2021-12-16";
###################################################################
my($cmd, $pathoList, $outdir, $version);
GetOptions(
	"cmd:s" => \$cmd,
	"pathoList:s" => \$pathoList,
	"outdir:s" => \$outdir,
	"version|v" => \$version,
);
if($version){
	print basename $0." $VERSION\n";
	exit(0);
}
&help unless($cmd && $pathoList);
$outdir ||= "./";
`mkdir -p $outdir` unless(-d $outdir);


###################################################################
my %k2name = (
	"Bacteria" => "bacteria",
	"Eukaryota:Fungi" => "fungi",
	"Viruses" => "viral",
	"Eukaryota:Parasite" => "protozoa",
	"Eukaryota:Protozoa" => "protozoa",
	"SpecialPathogens" => "LY",
);

###################################################################
open IN,"$pathoList" or die $!;
<IN>;
my %hash;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if($t[1] eq '-'){
	    $hash{$t[0]}{"NULL_$t[0]"} = 1;
	}else{
		$hash{$t[0]}{$t[1]} = 1;
	}
}
close IN;

open IN,"$cmd" or die $!;
open OU,">$outdir/BWA.sh" or die $!;
my(%pack, $path);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if(exists $hash{$t[0]}{$t[1]}){
		$path = $t[2];
		foreach my $i(split /,/,$t[3]){
			push @{$pack{$k2name{$t[0]}}}, $i;
		}
		print OU "$t[4]\n";
	}
}
close IN;
close OU;

#`cd $path`;
#foreach my $i(keys %pack){
#	my $line = join " ", @{$pack{$i}};
#	`cd $path && tar --use-compress-program=$PIGZ -cpf $path/$i.result.fq.gz $line`;
#
#}

sub help{
print "
		Usage: perl $0

		--cmd           list of command
		--pathoList     list of pathogeny
		--outdir        [./]

";
exit(0);
}
