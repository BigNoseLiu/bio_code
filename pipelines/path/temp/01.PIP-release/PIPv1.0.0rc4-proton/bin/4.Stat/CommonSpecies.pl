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
	my $file_nrow=1; #hy
	while(<FL>){
		chomp;
		my @s = split /\t/,$_;
		if($s[5] eq '-'){
			$hash->{$sample}->{$type}->{"$s[0]\t$s[3]\t-"}->{Read} = $s[7];
			$hash->{$sample}->{$type}->{"$s[0]\t$s[3]\t-"}->{RPM} = $s[10];
			$hash->{$sample}->{$type}->{"$s[0]\t$s[3]\t-"}->{Normalize} = $s[11];
			$hash2->{"$s[0]\t$s[3]\t-"} += 1;
		}else{
			$hash->{$sample}->{$type}->{"$s[0]\t$s[5]\t$s[6]"}->{Read} = $s[7];
			$hash->{$sample}->{$type}->{"$s[0]\t$s[5]\t$s[6]"}->{RPM} = $s[10];
			$hash->{$sample}->{$type}->{"$s[0]\t$s[5]\t$s[6]"}->{Normalize} = $s[11];
			$hash2->{"$s[0]\t$s[5]\t$s[6]"} += 1;
		}
		$file_nrow+=1; #hy
	}
	close FL;
	
	
	if($file_nrow eq 1){$hash->{$sample}->{$type}->{"none"}->{Read}=0; $hash->{$sample}->{$type}->{"none"}->{RPM}=0;$hash->{$sample}->{$type}->{"none"}->{Normalize}=0;}# hy    #Total或Pathogeny表只有title的样本，也输出到样本总表
	
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
				#if(exists $stat->{$i}->{$j}->{$k}->{Read}){ #hy
				if(exists $stat->{$i}->{$j}->{$k}){ #hy
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
	open OU3,">$out.Normalize_10M.txt" or die $!;
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
	
	# output 9 species -start hy 20220524
	my @speslistA=("Klebsiella pneumoniae","Acinetobacter baumannii","Pseudomonas aeruginosa","Staphylococcus aureus","Haemophilus influenzae","Streptococcus pneumoniae","Stenotrophomonas maltophilia","Escherichia coli","Aspergillus fumigatus");
	my %focus_species_here;
	foreach my $spej(@speslistA){
	$focus_species_here{$spej}++;
	}

	open IN,"$out.Reads.txt";
	open OUT,">$out.Reads_9species.txt";
	while(<IN>){
		chomp;
		if($_=~/^\s*$/){next;} 
		my @c=split /\t/,$_;
		if(($c[0] eq "Kingdom") or ($c[0] eq "-")){
			print OUT "$_\n";
		}elsif(exists $focus_species_here{$c[1]}){
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
	
	open IN,"$out.RPM.txt";
	open OUT,">$out.RPM_9species.txt";
	while(<IN>){
		chomp;
		if($_=~/^\s*$/){next;} 
		my @c=split /\t/,$_;
		if(($c[0] eq "Kingdom") or ($c[0] eq "-")){
			print OUT "$_\n";
		}elsif(exists $focus_species_here{$c[1]}){
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
	
	open IN,"$out.Normalize_10M.txt";
	open OUT,">$out.Normalize_10M_9species.txt";
	while(<IN>){
		chomp;
		if($_=~/^\s*$/){next;} 
		my @c=split /\t/,$_;
		if(($c[0] eq "Kingdom") or ($c[0] eq "-")){
			print OUT "$_\n";
		}elsif(exists $focus_species_here{$c[1]}){
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;

	
	`$PERL $TXT2EXCEL $out.Reads.txt $out.Reads.xlsx`;
	`$PERL $TXT2EXCEL $out.RPM.txt $out.RPM.xlsx`;
	`$PERL $TXT2EXCEL $out.Normalize_10M.txt $out.Normalize_10M.xlsx`;
	
	`$PERL $TXT2EXCEL $out.Reads_9species.txt $out.Reads_9species.xlsx`;
	`$PERL $TXT2EXCEL $out.RPM_9species.txt $out.RPM_9species.xlsx`;
	`$PERL $TXT2EXCEL $out.Normalize_10M_9species.txt $out.Normalize_10M_9species.xlsx`;
}

sub help{
print "
		Usage: perl $0

		--list        sample list
		--outdir      [./]

";
exit(0);
}
