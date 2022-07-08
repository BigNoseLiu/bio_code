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


my ( $help, $infile, $out_dir, $stat_flag);
$stat_flag = 1;
GetOptions(
	"help"=>\$help,
	"in=s"=>\$infile,
	"out=s"=>\$out_dir,
	"stat=s"=>\$stat_flag,
);

my %h_to_out_job = ();
`mkdir -p $out_dir`;

my $global_config = "";
my $global_config_flag = 0;
my $job_id = "NA";
my $last_job_ids = "NA";
my %h_jobs = ();
my %depend_jobs = ();
my $line;
while( $line = <STDIN> ){

	chomp $line;
	my @arr = split(/\t/,$line);
	#任务id
	if( $line =~ /^\s*#&J/ ){
		$global_config_flag = 0;
		$job_id = $arr[1];
		$last_job_ids = $arr[2];
	}
	#全局变量
	if( $line =~ /^\s*#&C/ ){
		$global_config_flag = 1;
	}

	if( $global_config_flag == 1 ){
		$global_config .= "$line\n";
	}
	elsif( $job_id ne "NA" ){
		next if( $job_id =~ /^AllStat/ && $stat_flag == 0 );
		$h_jobs{$job_id}{"shell"} .= $line."\n";
		foreach my $last_job_id( split(/,/,$last_job_ids) ){
			if( !defined( $h_jobs{$job_id}{"last_job"}{ $last_job_id } ) ){
				$h_jobs{$job_id}{"count"}++;
				$h_jobs{$job_id}{"last_job"}{ $last_job_id } = $h_jobs{$job_id}{"count"};
				$depend_jobs{ $last_job_id }++;
			}
		}
	}
}
foreach my $job_id( sort {$a cmp $b} keys(%h_jobs) ){
	if( !defined($depend_jobs{$job_id}) ){
		open OUT,">$out_dir/sh_$job_id.sh" or die $!;
		print OUT $global_config."\n".&add_depend_shell( $job_id,$job_id )."\n";
		close OUT;
	}
}

sub add_depend_shell{
	my ( $ori_id, $job_id ) = @_;
	my $out_str = "";
	if( !defined( $h_to_out_job{$ori_id}{$job_id} ) ){
		$out_str = $h_jobs{$job_id}{"shell"};
		$h_to_out_job{$ori_id}{$job_id}++;
		foreach my $last_job_id( sort { $h_jobs{$job_id}{"last_job"}{$a} <=> $h_jobs{$job_id}{"last_job"}{$b} } keys(%{$h_jobs{$job_id}{"last_job"}}) ){
			if( $last_job_id eq "NA" || defined($h_to_out_job{$ori_id}{$last_job_id}) ){
				next;
			}
			$out_str = &add_depend_shell($ori_id,$last_job_id)."\n$out_str";
			$h_to_out_job{$ori_id}{$last_job_id}++;
		}
	}
	return $out_str;
}
