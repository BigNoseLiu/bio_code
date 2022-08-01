#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use MongoDB;
use Encode;
use BSON;
#use BSON::Doc;
use BSON::Types ':all';


my $db_col_name = $ARGV[0];
my $task_id = $ARGV[1];


my $client = MongoDB::MongoClient->new(
   host => '10.10.9.35:31112',
   username => 'lintop',
   password => 'lintop.hx321.mongo'
);
$client->reconnect;
my $collection = $client->ns( $db_col_name );

my %h_stat = ();
my %h_sig= ();
my $line;
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	my ( $gene, $trans, $type, $sig, $count ) = @arr;
	if( $sig =~ /Pathogenic/ ){
		$sig = "Pathogenic";
	}
	elsif( $sig =~ /Likely_pathogenic/ ){
		$sig = "Likely_pathogenic";
	}
	elsif( $sig =~ /Benign/ ){
		$sig = "Benign";
	}
	elsif( $sig =~ /Likely_benign/ ){
		$sig = "Likely_benign";
	}
	elsif( $sig =~ /Uncertain_significance/ ){
		$sig = "Uncertain_significance";
	}
	else{
		next;
	}
	$h_stat{$gene}{$trans}{$type}{$sig} += $count;
	$h_sig{$sig}++;
}

foreach my $gene( sort {$a cmp $b} keys(%h_stat) ){
	foreach my $trans( sort {$a cmp $b} keys(%{$h_stat{$gene}}) ){
		my %h_record = ('gene'=>$gene,'transcript'=>$trans);
		my @arr_trans = ();
		foreach my $type( sort {$a cmp $b} keys(%{$h_stat{$gene}{$trans}}) ){
			foreach my $sig( sort {$a cmp $b} keys(%h_sig) ){
				my $count = 0;
				if(defined($h_stat{$gene}{$trans}{$type}{$sig})){
					$count = $h_stat{$gene}{$trans}{$type}{$sig};
				}
				$h_record{'effects'}{$type}{'clinvar_stat'}{$sig} = $count;
				$h_record{'effects'}{$type}{'ch_name'} = $type;
			}
		}
		$collection->insert_one( \%h_record );
	}
}



$client->disconnect;
