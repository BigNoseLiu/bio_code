#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use FindBin qw($Bin $Script);

my %h_fqs = ();
my $line;
while( $line = <STDIN> ){
	#202205231446_210901006_B_POOLINHG18_ZMCAFY4085_L01_R2.fq.gz
	chomp $line;
	if( $line =~ /(POOLINHG\d+_[^_]+)_[^_]+_R1/ ){
		my $sample=$1;
		if( !defined($h_fqs{ $sample }{'R1'}) ){
			$h_fqs{ $sample }{'R1'} = $line;
		}
		else{
			print STDERR "dup sample for $sample R2\n";
		}
	}
	elsif( $line =~ /(POOLINHG\d+_[^_]+)_[^_]+_R2/ ){
		my $sample=$1;
		if( !defined($h_fqs{ $sample }{'R2'}) ){
			$h_fqs{ $sample }{'R2'} = $line;
		}
		else{
			print STDERR "dup sample for $sample R2\n";
		}
	}
	else{
		print STDERR "illegal fomrat for $line\n";
	}
}





my $uid="0:0";
my $user_ids=`id`;
if( $user_ids =~ /^uid=(\d+)\(\S+\s+gid=(\d+)\(/ ){
	$uid="$1:$2";
}
else{
	print STDERR "Err:unable to get uid\n";
}

my $cwd =`pwd`;	chomp $cwd;
my $out_dir = "$cwd/out_data";
$out_dir = $ARGV[0] if(defined($ARGV[0]));
`mkdir -p $out_dir` if(!-d $out_dir);

foreach my $sample( sort {$a cmp $b} keys(%h_fqs) ){
	print "sh $Bin/HPV_16_18.v1.sh $uid $sample $out_dir/$sample ".$h_fqs{$sample}{'R1'}." ".$h_fqs{$sample}{'R2'}."\n";
}
