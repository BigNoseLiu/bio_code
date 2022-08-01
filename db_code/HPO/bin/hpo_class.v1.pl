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
my $chpo_file = $ARGV[1];
my %h_chpo = ();
my $line;
open IN,$chpo_file or die $!;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_chpo{$arr[1]}{'ch_name'} = decode_utf8($arr[4]) if(defined($arr[4]) && $arr[4] =~ /\S/);
	$h_chpo{$arr[1]}{'ch_def'} = decode_utf8($arr[5]) if(defined($arr[5]) && $arr[5] =~ /\S/);
}
close IN;


my $client = MongoDB::MongoClient->new(
   host => '1.12.236.215:27017',
   username => 'lintop',
   password => 'lintop.hx321.mongo'
);
$client->connect;
#my $client = MongoDB->connect( '10.10.9.22:27017' );


my $collection = $client->ns( $db_col_name );

my %h_record= ();
my @arr_syn = ();
while( $line = <STDIN> ){
	chomp $line;
	if( $line =~ /^\s*\[Term\]\s*$/ ){
		if( defined($h_record{'id'}) ){
			if(scalar@arr_syn > 0){
				$h_record{'synonym'} = \@arr_syn;
			}
			$collection->insert_one( \%h_record );
		}
		%h_record= ();
		@arr_syn = ();
	}
	elsif($line =~ /^id: HP:(\d+)\s*$/){
		my $id = $1;
		$h_record{'id'} = "HP_$id";
		$h_record{'hpo_id'} = "HP:$id";
		if(defined($h_chpo{"HP:$id"}{'ch_name'})){
			$h_record{'ch_name'} = $h_chpo{"HP:$id"}{'ch_name'};
		}
		if(defined($h_chpo{"HP:$id"}{'ch_def'})){
			$h_record{'ch_def'} = $h_chpo{"HP:$id"}{'ch_def'};
		}
	}
	elsif($line =~ /^name: (.+\S)\s*$/){
		$h_record{'en_name'} = $1;
	}
	elsif($line =~ /^def: "([^"]+)"/){
		$h_record{'en_def'} = $1;
	}
	elsif($line =~ /^alt_id: HP:(\d+)\s*$/){
		$h_record{'alt_id'}{"HP_$1"} = "alt_id";
	}
	elsif($line =~ /^is_a: HP:(\d+)\s/){
		$h_record{'parent'}{"HP_$1"} = "is_a";
	}
	elsif($line =~ /^synonym: "([^"]+)"/){
		my $synonym = $1;
		my %h_syn = ();
		$h_syn{'en_name'} = $synonym;
		push @arr_syn, \%h_syn;
	}
	elsif($line =~ /^xref: UMLS:(\S+)\s*$/){
		$h_record{'UMLS'}{$1} = "xref";
	}
	elsif($line =~ /^xref: MSH:(\S+)\s*$/){
		$h_record{'MSH'}{$1} = "xref";
	}
	elsif($line =~ /^xref: SNOMEDCT_US:(\S+)\s*$/){
		$h_record{'SNOMEDCT_US'}{$1} = "xref";
	}
	my @arr = split(/\t/,$line);
}
if( defined($h_record{'id'}) ){
	$collection->insert_one( \%h_record );
}

$client->disconnect;
$client2->disconnect;
