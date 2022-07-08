#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor;

sub usageWithColor {
	print color "green";#change the text color
	print <<USAGE;
Description:
        Used to normalize scores for employees.
Usage:  
        perl $0  -in in_file -out out_file
        Options:
	     -help : reveal help info.
	     -in   <str>  : in file path.
		File format : 
			field 1: name of employees.
			field 2: basic performance payment of employees.
			field 3: performance judged by each group leader.
				Note: this stand for group score judged by yanzhixiang when comes to the group leader's row, and the score of group leader is defaulted as 90.
			field 4: group id of employees 
				Example: 1.2 means the No.2 member of No.1 group.
				Note: group leader must be X.0 which X stand for group ID.
	     -out  <str>  : output file path
			field 3: this is the final score reported.
			other fields are the same with those of in_file.
Example:
       perl $0 -in score.txt  -out  output_file
Author & Contact:
	Mingming Liu
        liumingming\@genomics.cn
Last updated:
        2013-12-11
USAGE
	print color "reset";#change back the text color
}

my ($in_file,$help,$out,$total_pay);
GetOptions(
	"help"=>\$help,
	"out=s"=>\$out,
	"in=s"=>\$in_file,
	#"totalPay=i"=>\$total_pay,
);
if ( defined $help || !defined( $out ) || !defined( $in_file ) ) {
	&usageWithColor();
	exit 0;
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}
open IN,"$in_file" or die "Fail opening $in_file\n";
my $line;
my $line_count = 0;
my (%h_name,%h_basic_perf,%h_perf,%h_group);
my (%h_group2payTotal,%h_perf_real,%h_group_real,%h_group2Score,%h_group2perfTotal);
while( $line = <IN> ){
	next if( $line !~ /\d/ );#skip the title
	$line_count++;
	chomp $line;
	$line =~ s/^\s+(\S)/$1/;
	my @arr = split(/\s+/,$line);

	#format checking
	#last 3 colomns must be numberic
	if( $arr[1] !~ /^[\d\.]+$/ || $arr[2] !~ /^[\d\.]+$/ || $arr[3] !~ /^[\d\.]+$/  ){
		&err_print("Error: last 3 colomns must be numberic.\nline $line_count :\t$line\n");
		exit 0;
	}
	#erase the space at the head and tail of each field
	for(my $t_i=0; $t_i < scalar@arr; $t_i++){
		$arr[$t_i] =~ s/^\s*(\S)/$1/;
		$arr[$t_i] =~ s/(\S)\s*$/$1/;
	}

	$h_name{$line_count}=$arr[0];		#name of a person
	$h_perf{$line_count}=$arr[1];		#performance judged by group leader
	$h_group{$line_count}=$arr[2];		#group id
	$h_basic_perf{$line_count}=$arr[3];	#basic performance pay

	$total_pay += $arr[3];			#get total perf pay
	#get the group info
	my ( $group_id,$mem_id );
	if( $arr[2] =~ /^(\d+)\.(\d+)$/ ){
		( $group_id,$mem_id ) = ($1,$2);
		$h_group_real{$line_count} = $group_id;
	}
	else{
		&err_print("Error: group_id wrong formated.\nline $line_count :\t$line\n");
		exit 0;
	}
	my $mem_score = $arr[1];
	if( $mem_id == 0 ){
		$mem_score = 90;
		$h_group2Score{$group_id} = $arr[1];     	#group score
	}

	$h_perf_real{$line_count} = $mem_score;
	$h_group2payTotal{$group_id} += $arr[3];		#basic performance pay of each group
	$h_group2perfTotal{$group_id} += ($arr[3]*$mem_score);	#denominator when calculate the ratio of each employ in their group

}
close IN;

my $total_groupPerf = 0;
my %h_group2perfPay;
foreach my $group_name( keys(%h_group2payTotal) ){
	if( defined( $h_group2payTotal{$group_name} ) && defined( $h_group2Score{$group_name} ) ){
		$total_groupPerf += (  $h_group2payTotal{$group_name} * $h_group2Score{$group_name}   );
	}
}

#calculate the payment of each group
foreach my $group_name( keys(%h_group2payTotal) ){
	if( defined( $h_group2payTotal{$group_name} ) && defined( $h_group2Score{$group_name} ) ){
		$h_group2perfPay{$group_name} = (($h_group2payTotal{$group_name} * $h_group2Score{$group_name})/$total_groupPerf)*$total_pay;
	}
}

open OUT,'>',"$out" or die "Fail opening $out\n";
print OUT "#Name\tFinal_score\tPerf_score\tGroup\tBasic_perf_pay\n";

my $actual_pay = 0;
my $max_final_perf = 0;
my %h_final_perf;
foreach my $count( sort {$a <=> $b} keys( %h_name ) ){
	my $final_perf = sprintf("%.2f",($h_perf_real{$count}/$h_group2perfTotal{ $h_group_real{$count} })*($h_group2perfPay{ $h_group_real{$count} }) );
	$actual_pay += $h_basic_perf{$count} * $final_perf;
	$final_perf = $final_perf*100;
	if( $final_perf > $max_final_perf ){
		$max_final_perf = $final_perf;
	}
	$h_final_perf{ $count } = $final_perf;
}

#normalize the score
foreach my $count( sort {$a <=> $b} keys( %h_name ) ){
	my $final_perf = int( ($h_final_perf{ $count }*100)/$max_final_perf );
	print OUT $h_name{$count}."\t".$final_perf."\t".$h_perf{$count}."\t".$h_group{$count}."\t".$h_basic_perf{$count}."\n";
}

close OUT;
