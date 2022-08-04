#!/usr/bin/env perl

#===============================================================================
#   FileName: QC.stat.pl
#   Author  : qiuwancen
#   E-Mail  : 972538446@qq.com
#   Version : 1.0
#   Date    : 2022-02-28
#===============================================================================

use strict;
use File::Basename;
use Cwd qw(abs_path);

die "perl $0 [QC.sample.list] [out]\n" unless(@ARGV==2);

#my $path = abs_path $ARGV[0];
my $path = $ARGV[0];
my $dir = dirname $path;
$dir = abs_path "$dir/..";

my %trans = (
	"total reads" => "Reads",
	"total bases" => "Bases",
	"Q20 bases" => "Bases Q20",
	"Q30 bases" => "Bases Q30",
	"reads failed due to low quality" => "Reads with Low Quality",
	"reads failed due to too many N" => "Reads with N",
	"reads failed due to too short" => "Reads filted due to too short",
	"reads failed due to low complexity" => "Read with Low Complexity",
	"reads with adapter trimmed" => "Reads with Adapter",
	"bases trimmed due to adapters" => "Bases with Adapter",
	"Duplication rate (may be overestimated since this is SE data)" => "Duplication Rate",
);

my @infos = ("Raw Reads", "Raw Bases", "Raw Bases Q20", "Raw Bases Q30", "Clean Reads", "Clean Bases", "Clean Bases Q20", "Clean Bases Q30", "Duplication Rate", "Reads with Adapter", "Bases with Adapter", "Reads with N", "Reads with Low Quality", "Read with Low Complexity", "Reads filted due to too short", "Host Reads", "Host Rate", "Internal Control Reads", "Internal Control Rate", "Internal Control Gene18 Average Depth", "Internal Control Gene18 Coverage", "Internal Control Gene43 Average Depth", "Internal Control Gene43 Coverage", "Classified_Reads_Number", "Classified_Rate", "Unclassified_Reads_Number", "Unclassified_Rate", "Viruses_Reads_Number", "Viruses_Rate", "Bacteria_Reads_Number", "Bacteria_Rate", "Fungi_Reads_Number", "Fungi_Rate", "Protozoa_Reads_Number", "Protozoa_Rate", "Metazoa_Parasite_Reads_Number", "Metazoa_Parasite_Rate");
open IN,"$path" or die $!;
my %hash;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if(-d "$dir/01.QC/01.Fastp"){
		open FP,"$dir/01.QC/01.Fastp/$t[0]/$t[1]/$t[0].log" or die $!;
		my($tmp, $flag);
		while(<FP>){
			chomp;
			next unless(/:/);
			next if(/report|WARNING|fastp|reads passed filter:/);
			if(/before/){
				$flag = "Raw";
			}elsif(/after/){
				$flag = "Clean";
			}
			if(/^Read1/){
				$tmp = "Reads";
			}elsif(/^Read2/){
				$tmp = "Reads";
			}elsif(/^Filtering/){
				$tmp = "";
				$flag = "";
			}
			if(/: /){
				my @s = split /: /,$_;
				if($tmp ne ""){
					if($s[1] =~ /[\(\)]/){
						$s[1] =~ s/.*\(|%\)//g;
						$s[1] = sprintf("%.2f%", $s[1]);
					}
					$hash{$t[0]}{$t[1]}{"$flag $trans{$s[0]}"} = $s[1];
				}else{
					if($s[1] =~ /%/){
						$s[1] =~ s/%.*//g;
						$s[1] = sprintf("%.2f%", $s[1]);
					}
					$hash{$t[0]}{$t[1]}{"$trans{$s[0]}"} = $s[1];
				}
			}
		}
		close FP;
	}
	if(-d "$dir/01.QC/02.RemoveHost"){
		open FP,"$dir/01.QC/02.RemoveHost/$t[0]/$t[1]/$t[0].Homo.log" or die $!;
		my($hostReads, $hostRate);
		while(<FP>){
			chomp;
			if(/ mapped /){
				($hostReads, $hostRate) = ($1, $2) if(/(\d+) \+ \S+ mapped \((\S+) : /);
				$hash{$t[0]}{$t[1]}{"Host Reads"} = $hostReads;
				last;
			}
		}
		close FP;
		if(-e "$dir/01.QC/02.RemoveHost/$t[0]/$t[1]/$t[0].delIT.log"){
			open FP,"$dir/01.QC/02.RemoveHost/$t[0]/$t[1]/$t[0].delIT.log" or die $!;
			my $ITReads;
			while(<FP>){
				chomp;
				if(/ mapped /){
					$ITReads = $1 if(/(\d+) \+ \S+ mapped \((\S+) : /);
					$hash{$t[0]}{$t[1]}{"Internal Control Reads"} = $ITReads;
					$hash{$t[0]}{$t[1]}{"Internal Control Rate"} = sprintf("%.2f%", $ITReads/$hash{$t[0]}{$t[1]}{"Raw Reads"}*100);
					last;
				}
			}
			close FP;
		}
		if(-e "$dir/01.QC/02.RemoveHost/$t[0]/$t[1]/region.tsv.gz"){
			open FP,"gzip -dc $dir/01.QC/02.RemoveHost/$t[0]/$t[1]/region.tsv.gz |" or die $!;
			while(<FP>){
				chomp;
				my @a = split /\t/,$_;
				if($a[1] == 2297){
					$hash{$t[0]}{$t[1]}{"Internal Control Gene18 Average Depth"} = $a[3];
					$hash{$t[0]}{$t[1]}{"Internal Control Gene18 Coverage"} = $a[5];
				}
				if($a[1] == 4238){
					$hash{$t[0]}{$t[1]}{"Internal Control Gene43 Average Depth"} = $a[3];
					$hash{$t[0]}{$t[1]}{"Internal Control Gene43 Coverage"} = $a[5];
				}
			}
			close FP;
		}
	}
	if(-d "$dir/02.Annotation/01.Taxonomy"){
		open FP,"$dir/02.Annotation/01.Taxonomy/$t[0]/$t[1]/Taxonomy_Summary.txt" or die $!;
		<FP>;
		while(<FP>){
			chomp;
			my @a = split /\t/,$_;
			$hash{$t[0]}{$t[1]}{Classified_Reads_Number} = $a[4] - $a[12];
			$hash{$t[0]}{$t[1]}{Classified_Rate} = sprintf("%.3f%", $hash{$t[0]}{$t[1]}{Classified_Reads_Number} / $hash{$t[0]}{$t[1]}{"Clean Reads"} * 100);#"$a[5]%";
			$hash{$t[0]}{$t[1]}{"Host Reads"} += $a[12];
			$hash{$t[0]}{$t[1]}{"Host Rate"} = sprintf("%.3f%", $hash{$t[0]}{$t[1]}{"Host Reads"} / $hash{$t[0]}{$t[1]}{"Clean Reads"} * 100);
			$hash{$t[0]}{$t[1]}{Unclassified_Reads_Number} = $a[2];
			$hash{$t[0]}{$t[1]}{Unclassified_Rate} = sprintf("%.3f%", $hash{$t[0]}{$t[1]}{Unclassified_Reads_Number} / $hash{$t[0]}{$t[1]}{"Clean Reads"} * 100);
			$hash{$t[0]}{$t[1]}{Viruses_Reads_Number} = $a[6];
			$hash{$t[0]}{$t[1]}{Viruses_Rate} = sprintf("%.3f%", $hash{$t[0]}{$t[1]}{Viruses_Reads_Number} / $hash{$t[0]}{$t[1]}{"Clean Reads"} * 100);
			$hash{$t[0]}{$t[1]}{Bacteria_Reads_Number} = $a[7];
			$hash{$t[0]}{$t[1]}{Bacteria_Rate} = sprintf("%.3f%", $hash{$t[0]}{$t[1]}{Bacteria_Reads_Number} / $hash{$t[0]}{$t[1]}{"Clean Reads"} * 100);
			$hash{$t[0]}{$t[1]}{Fungi_Reads_Number} = $a[9];
			$hash{$t[0]}{$t[1]}{Fungi_Rate} = sprintf("%.3f%", $hash{$t[0]}{$t[1]}{Fungi_Reads_Number} / $hash{$t[0]}{$t[1]}{"Clean Reads"} * 100);
			$hash{$t[0]}{$t[1]}{Protozoa_Reads_Number} = $a[10];
			$hash{$t[0]}{$t[1]}{Protozoa_Rate} = sprintf("%.3f%", $hash{$t[0]}{$t[1]}{Protozoa_Reads_Number} / $hash{$t[0]}{$t[1]}{"Clean Reads"} * 100);
			$hash{$t[0]}{$t[1]}{Metazoa_Parasite_Reads_Number} = $a[11];
			$hash{$t[0]}{$t[1]}{Metazoa_Parasite_Rate} = sprintf("%.3f%", $hash{$t[0]}{$t[1]}{Metazoa_Parasite_Reads_Number} / $hash{$t[0]}{$t[1]}{"Clean Reads"} * 100);
		}
		close FP;
	}
}
close IN;

open OU,">$ARGV[1]" or die $!;
my($head1, $head2) = ("Sample", "Type");
my %anno;
foreach my $i(sort keys %hash){
	$head1 .= "\t$i";
	foreach my $j(sort keys %{$hash{$i}}){
		$head2 .= "\t$j";
		foreach my $k(@infos){
			$anno{$k} .= "\t$hash{$i}{$j}{$k}";
		}
	}
}
print OU "$head1\n$head2\n";
foreach my $i(@infos){
	print OU "$i$anno{$i}\n";
}
close OU;