#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use FindBin qw($Bin $Script);

my %h_fqs = ();
my $line;
while( $line = <STDIN> ){
	#202207191826_220101013_B_PE75_DJY0719_46_L01_R2.fq.gz
	chomp $line;
	if( $line =~ /(DJY\d+_[^_]+)_[^_]+_R1/ ){
		my $sample=$1;
		if( !defined($h_fqs{ $sample }{'R1'}) ){
			$h_fqs{ $sample }{'R1'} = $line;
		}
		else{
			print STDERR "dup sample for $sample R2\n";
		}
	}
	elsif( $line =~ /(DJY\d+_[^_]+)_[^_]+_R2/ ){
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

my $cwd =`pwd`;	chomp $cwd;
my $out_dir = "$cwd/out_data";
$out_dir = $ARGV[0] if(defined($ARGV[0]));
my $uid="0:0";
my $user_ids=`id`;
if( $user_ids =~ /^uid=(\d+)\(\S+\s+gid=(\d+)\(/ ){
	$uid="$1:$2";
}
else{
	print STDERR "Err:unable to get uid\n";
}

`mkdir -p $out_dir` if(!-d $out_dir);
foreach my $sample( sort {$a cmp $b} keys(%h_fqs) ){
	print "sh $Bin/DJY_UMI_PE.v1_02.sh  $uid $sample $out_dir 75T ".$h_fqs{$sample}{'R1'}." 7M1S67T ".$h_fqs{$sample}{'R2'}."\n";
}
