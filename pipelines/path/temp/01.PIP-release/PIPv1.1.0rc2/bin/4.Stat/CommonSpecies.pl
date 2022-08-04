#!/usr/bin/env perl
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename;
use Cwd qw(abs_path);
use lib "$Bin/../../lib/perl5";
use MyModule::MyLog qw(printInfo2ERR);
use MyModule::GlobalVar qw($TXT2EXCEL $PERL);

my($list, $outdir);
GetOptions(
	"list:s" => \$list,
	"outdir:s" => \$outdir,
);
&help unless($list);
$outdir ||= "./";
$outdir = abs_path $outdir;
`mkdir -p $outdir` unless(-e $outdir);

open IN,"$list" or die $!;
my(%total, %patho, %spe);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	&getCommon($t[-2], $t[0], $t[1], \%total, \%spe);
	&getCommon($t[-1], $t[0], $t[1], \%patho, \%spe);
}
close IN;

&Stat(\%total, "$outdir/Total_CommonSpecies");
&Stat(\%patho, "$outdir/Pathogeny_CommonSpecies");

sub getCommon{
	my($report, $sample, $type, $hash, $hash2) = @_;
	open FL,"$report" or die $!;
	<FL>;
	while(<FL>){
		chomp;
		my @s = split /\t/,$_;
		if($s[5] eq '-'){
			$hash->{$sample}->{$type}->{"$s[0]\t$s[7]\t-"}->{Read} = $s[13];
			$hash->{$sample}->{$type}->{"$s[0]\t$s[7]\t-"}->{RPM} = $s[18];
			$hash->{$sample}->{$type}->{"$s[0]\t$s[7]\t-"}->{Normalize} = $s[14];
			$hash2->{"$s[0]\t$s[7]\t-"} += 1;
		}else{
			$hash->{$sample}->{$type}->{"$s[0]\t$s[11]\t$s[12]"}->{Read} = $s[13];
			$hash->{$sample}->{$type}->{"$s[0]\t$s[11]\t$s[12]"}->{RPM} = $s[18];
			$hash->{$sample}->{$type}->{"$s[0]\t$s[11]\t$s[12]"}->{Normalize} = $s[14];
			$hash2->{"$s[0]\t$s[11]\t$s[12]"} += 1;
		}
	}
	close FL;
}

sub Stat{
	my($stat, $out) = @_;
	my(%sum, %count, %reads, %RPM, %Norma);
	my @Samples = sort keys %$stat;
	my $types;
	my $sampleNum;
	foreach my $i(@Samples){
		foreach my $j(sort keys %{$stat->{$i}}){
			$sampleNum += 1;
			$types .= "\t$j";
			foreach my $k(keys %spe){
				if(exists $stat->{$i}->{$j}->{$k}->{Read}){
					$sum{$k} += $stat->{$i}->{$j}->{$k}->{Read};
					$count{$k} += 1;
					$reads{$k} .= "\t$stat->{$i}->{$j}->{$k}->{Read}";
					$RPM{$k} .= "\t$stat->{$i}->{$j}->{$k}->{RPM}";
					$Norma{$k} .= "\t$stat->{$i}->{$j}->{$k}->{Normalize}";
				}else{
					$reads{$k} .= "\t0";
					$RPM{$k} .= "\t0";
					$Norma{$k} .= "\t0";
				}
			}
		}
	}
	open OU1,">$out.Reads.txt" or die $!;
	open OU2,">$out.RPM.txt" or die $!;
	open OU3,">$out.Normalize.txt" or die $!;
	print OU1 "Kingdom\tScientificName\tChineseName\tComment\t".(join "\t", @Samples)."\n";
	print OU1 "-\t-\tType\t-$types\n";
	print OU2 "Kingdom\tScientificName\tChineseName\tComment\t".(join "\t", @Samples)."\n";
	print OU2 "-\tType\t-$types\n";
	print OU3 "Kingdom\tScientificName\tChineseName\tComment\t".(join "\t", @Samples)."\n";
	print OU3 "-\t-\tType\t-$types\n";
	foreach my $keys(sort {$sum{$b} <=> $sum{$a}} keys %sum){
		print OU1 "$keys\t$count{$keys}/$sampleNum$reads{$keys}\n";
		print OU2 "$keys\t$count{$keys}/$sampleNum$RPM{$keys}\n";
		print OU3 "$keys\t$count{$keys}/$sampleNum$Norma{$keys}\n";
	}
	close OU1;
	close OU2;
	close OU3;
	`$PERL $TXT2EXCEL $out.Reads.txt $out.Reads.xlsx`;
	`$PERL $TXT2EXCEL $out.RPM.txt $out.RPM.xlsx`;
	`$PERL $TXT2EXCEL $out.Normalize.txt $out.Normalize.xlsx`;
}

sub help{
print "
		Usage: perl $0

		--list        sample list
		--outdir      [./]

";
exit(0);
}
