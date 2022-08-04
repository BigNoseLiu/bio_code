#!/usr/bin/env perl

#===============================================================================
#   FileName: m82annoTaxid.pl
#   Author  : 972538446@qq.com
#   Version : 1.0
#   Date    : 2021-06-02
#   Description: The program to annotate m8.
#===============================================================================

use warnings;
use strict;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use Getopt::Long;
use Pod::Usage;
use lib "$Bin/../lib";
use MyModule::GlobalVar qw($BIN_PATH $ID2TAX $SEQID2SPEID);
use MyModule::MyLog qw(printInfo2ERR);

my(
	$inM8, 
	$out, 
	$help
);

GetOptions(
	"inM8:s" => \$inM8,
	"out:s" => \$out,
	"help|h" => \$help,
);

if($help || ! $inM8)
{
	pod2usage(
		-exitval => 1, 
		-verbose => 1
	);
}

#===============================================================================
#   Set the default values of options
#===============================================================================
$out ||= "./out.txt";
$out = abs_path $out;

#===============================================================================
#   Check file
#===============================================================================
if(!-e $ID2TAX)
{
	die "Error: $ID2TAX does not exist.\n";
}
else
{
	printInfo2ERR("Table used to convert taxid into taxonomy: $ID2TAX");
}
if(!-e $SEQID2SPEID)
{
	die "Error: $SEQID2SPEID does not exist.\n";
}
else
{
	printInfo2ERR("Table used to convert sequence taxid into species taxid: $SEQID2SPEID");
}

#===============================================================================
#   Read comment file
#===============================================================================
printInfo2ERR("Reading comment file");
open IN,"$ID2TAX" or die $!;
<IN>;
my %taxid2tax;
while(<IN>)
{
	chomp;
	my @t = split /\t/,$_;
	$t[0] =~ s/ //g;
	$taxid2tax{$t[0]} = "$t[1]\t$t[-1]";
}
close IN;

open IN,"$SEQID2SPEID" or die $!;
my %seqid2speid;
while(<IN>)
{
	next if(/NULL/);
	chomp;
	my @t = split /\t/,$_;
	$seqid2speid{$t[0]} = $t[1];
}
close IN;


#===============================================================================
#   Read input file and output comment results
#===============================================================================
printInfo2ERR("Input file: $inM8");
printInfo2ERR("Output file: $out");
printInfo2ERR("Start commenting");
open IN,"$inM8" or die $!;
open OU,">$out" or die $!;
print OU "qseqid\tsacc\tpident\tlength\tevalue\tbitscore\tqcovhsp\tspe taxid\tkingdom\tspecies\n";
while(<IN>)
{
	chomp;
	my @t = split /\t/,$_;
	my $queryTaxID = (split /\|/,$t[0])[0];
	my $targetTaxID;
	if($t[1] =~ /:/)
	{
		$targetTaxID = (split /:/,$t[1])[0];
	}
	else
	{
		$targetTaxID = (split /\|/,$t[1])[0];
	}
	my($targetSpeID, $targetSpe) = ('NULL', 'NULL');
	if(exists $seqid2speid{$targetTaxID})
	{
		$targetSpeID = $seqid2speid{$targetTaxID};
	}
	if(exists $taxid2tax{$targetSpeID})
	{
		$targetSpe = $taxid2tax{$targetSpeID};
	}
	print OU "$_\t$targetSpeID\t$targetSpe\n";
}
close IN;
close OU;
printInfo2ERR("Finish commenting");

__END__


=pod

=head1 NAME

m82annoTaxid.pl - Annotate m8.

=head1 VERSION

V1.0

=head1 SYNOPSIS

perl m82annoTaxid.pl [options]

=head1 OPTIONS

=over 8

=item B<--inM8>

Input m8 file.

=item B<--out> <Path>

output path. [./]

=item B<--help|h>

Print this information.

=back

=cut
