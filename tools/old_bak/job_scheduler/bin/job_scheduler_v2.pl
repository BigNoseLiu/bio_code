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


#global settings
my $config_sperator = ":";

sub usageWithColor {
	print color "green";#change the text color
	print <<USAGE;
Description:
        Used to schedule all the steps by runing them in proper order to finsh a whole pipeline.
Usage:  
        perl $0 [options]
        Options:
	     -help : reveal help info
	     -schedule  <str>  : schedule file which illustrate all the steps to be done in order of the pipeline [required]
	     	Format:
			#Actived: true/false
			1 $config_sperator 1.1 1.2
			qsub=qsub -cwd -l vf=2g  step1.sh

			1.1 $config_sperator
			sh=sh step1.1.sh

			1.2 $config_sperator
			qsub=qsub -cwd -l vf=1g step1.2.sh

			2.1 $config_sperator
			sh=sh step2.1.sh

			2 $config_sperator 2.1
			qsub=qsub -cwd -l vf=1g step2.sh

	     -log  <str>  : log file to check out the current step of the whole pipeline [required]
			    only the lines with '^#' will be detected to check the done jobs
	     		    can be edited to re-run the jobs with bugs or unfinished
	     	Format:
			#1.1 done
			2014-1-10 16:20:10~~~2014-2-10 16:22:10
			sh step1.1.sh

			#1.2 done
			2014-1-10 16:20:10~~~2014-2-10 16:22:10
			qsub -cwd -l vf=1g step1.2.sh

			#1 done
			2014-1-10 16:20:10~~~2014-2-10 16:22:10
			qsub -cwd -l vf=2g step1.sh
        Example:
             perl $0 -s test.schedule -l log
Author & Contact:
	Mingming Liu
        mingming_liu\@shengtinggroup.com
Last updated:
        2014-07-09
USAGE
	print color "reset";#change back the text color
}

my ($help,$schedule_file,$log_file);
my $duration = 500;
GetOptions(
	"help"=>\$help,
	"schedule=s"=>\$schedule_file,
	"log=s"=>\$log_file,
	"duration=s"=>\$duration,
);
if (defined $help || (!defined $schedule_file) || (!defined $log_file) ) {
	&usageWithColor();
	#&test();
	exit 0;
}
open CONFIG,"$Bin/qsub.config" or die $!;
my $qsub_prefix = <CONFIG>;
chomp $qsub_prefix;
close CONFIG;

my $line;
my %h_done;
my $log_str = "";
if( -e $log_file ){
	open LOG,"$log_file" or die "Fail opening $log_file\n";
	while( $line = <LOG> ){
		chomp $line;
		next if( $line =~ /^\s*$/ );
		if( $line =~ /^#(\S+)\s+done\s*$/ ){
			$h_done{$1} = "$line\n";#done information
			$line = <LOG>;
			$h_done{$1} .= "$line";#time stamp
			$line = <LOG>;
			$h_done{$1} .= "$line";#cmd details
		}
	}
	close LOG;
}

open LOG,">>$log_file" or die "Fail opening $log_file\n";
$| = 0;
my %h_qsub_running;
my $done_flag = 0;
while($done_flag == 0){

	#indicating if the pipeline is done
	$done_flag = 1;

	#check the qsub jobs running
	my $qstat_jobs = `qstat|sed '1,2d'|awk '{print \$1}'`;#system('qstat|sed \'1,2d\'|awk '{print \$1}'');
	my %h_tmp = ();
	foreach my $job_id( split(/\n/,$qstat_jobs) ){
		$h_tmp{$job_id}++;
	}

	#mark the jobs done already
	foreach my $target_id( keys( %h_qsub_running ) ){
		#jod done
		if( !defined( $h_tmp{ $h_qsub_running{$target_id}{"jobID"} } ) ){
			$h_done{$target_id}++;
			print LOG "#$target_id done\n";
			print LOG $h_qsub_running{$target_id}{"startTime"}."~~~".&getTime()."\n";
			print LOG $h_qsub_running{$target_id}{"cmd"}."\n";
			delete( $h_qsub_running{$target_id} );
		}
	}

	#my $actived = "false";
	open SCH,"$schedule_file" or die "Fail opening $schedule_file\n";
	my $cmd_line;
	#active information [Must be the first line]
	#$line = <SCH>;
	#if( $line =~ /^#Actived:/ ){
	#	if( $line =~ /true/i ){
			#$actived = "true";
			while( $line = <SCH> ){
				chomp $line;
				next if( $line =~ /^#/ );#annotaions lines
				next if( $line =~ /^\s*$/ );#blank lines

				#get the target ID & prerequisites of one step
				my ($tmp_id,$tmp_prerequisites);
				if( $line =~ /^\s*(\S+)\s*$config_sperator(.*)$/ ){
					$tmp_id = $1;
					$tmp_prerequisites = $2;
				}
				else{
					print STDERR "Wrong-Format-Schedule:prerequisites not correct\n$line\n";
					exit;
				}
				$cmd_line = <SCH>;#cmd line of this job

				#check if the job done already
				if( $h_done{$tmp_id} ){
					next;
				}
				else{
					$done_flag = 0;
				}
				#check if the job is already running
				next if( defined( $h_qsub_running{$tmp_id} ) );

				#check if the prerequisites done already
				my $prerequisitesDone_flag = 1;
				if( defined($tmp_prerequisites) && $tmp_prerequisites !~ /^\s*$/ ){
					#erase the spaces flanking the line
					$tmp_prerequisites =~ s/^\s*(\S)/$1/;
					$tmp_prerequisites =~ s/(\S)\s*$/$1/;
					foreach my $id( split(/\s+/,$tmp_prerequisites) ){
						if( !defined(  $h_done{$id} ) ){
							$prerequisitesDone_flag = 0;
							last;
						}
					}
				}
				#run the job if all the prerequisites are done
				if( $prerequisitesDone_flag == 1 ){
					&run_job( $tmp_id,$cmd_line );
				}
			}
	#	}
	#}
	#else{	#wrong formated lines
	#	print STDERR "Wrong-Format-Schedule:First line not active information\n";
	#	print STDERR "Please put "."#Actived:true"." in the first line\n";
	#	exit;
	#}

	close SCH;

	#check the status of the jobs every 8 minutes
	sleep($duration);
}


#submit the job & get the information of the job [including jobID,startTime,cmd used to submit]
sub run_job{
	my ($tmp_id,$tmp_cmdLine) = @_;
	#my $qsub_cmd = "qsub -cwd -q ngb.q -P ngb_db -l vf=";
	if( $tmp_cmdLine =~ /^qsub=(\S.*)$/ ){	#submit the job intended to be run on HPC
		my $qsub_cmd = "$qsub_prefix $1";
		#$qsub_cmd .= $1;
		#qsub the job
		my $tmp_str = `$qsub_cmd`;
		#get the job ID being qsubed
		if( defined($tmp_str) && $tmp_str =~ /Your job\s+(\d+)\s/  ){
			#get the qsub ID of the job
			$h_qsub_running{$tmp_id}{"jobID"} = $1;
			$h_qsub_running{$tmp_id}{"startTime"} = &getTime();
			$h_qsub_running{$tmp_id}{"cmd"} = $qsub_cmd;
		}
		else{
			#qsub failed
			print STDERR "Fail at qsub : $qsub_cmd\n";
			exit;
		}
	}
	elsif( $tmp_cmdLine =~ /^sh=(\S.*)$/ ){	#run the cmd locally
		my $sh_cmd = $1;
		my $tmp_startTime = &getTime();
		`$sh_cmd`;
		$h_done{$tmp_id}++;
		print LOG "#$tmp_id done\n";
		print LOG $tmp_startTime."~~~".&getTime()."\n";
		print LOG "$sh_cmd\n";
	}
	else{
		#cmd line forat is wrong
		print STDERR "Wrong-Format-Schedule:cmd line not correct\n$line\n";
		exit;
	}
}

close LOG;


#get the time in a proper format
sub getTime{
	my ($week,$mon,$day,$ht,$year) = split(" ",localtime(time()) );
	return join("-",$year,$mon,$day)." $ht";
}
sub test{
	#check if the interface system will block the perl script
	#result: yes, it will
	print "1-".&getTime()."\n";
	system("sleep 10");
	print "2-".&getTime()."\n";
	`sleep 10`;
	print "3-".&getTime()."\n";
}
