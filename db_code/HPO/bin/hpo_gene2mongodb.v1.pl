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

my %h_hpo2gene = ();
my $line;

my $phenotype_to_genes = $ARGV[1];
open IN,$phenotype_to_genes or die $!;
while( $line = <IN> ){
	chomp $line;
	next if( $line =~ /^#/ );
	my @arr = split(/\t/,$line);
	$h_hpo2gene{$arr[0]}{$arr[3]}{'type'} = 'parent';
	$h_hpo2gene{$arr[0]}{$arr[3]}{'ref'}{$arr[6]}++;
}
close IN;

my $genes_to_phenotype = $ARGV[2];
open IN,$genes_to_phenotype or die $!;
while( $line = <IN> ){
	chomp $line;
	next if( $line =~ /^#/ );
	my @arr = split(/\t/,$line);
	$h_hpo2gene{$arr[2]}{$arr[1]}{'type'} = 'precise';
	$h_hpo2gene{$arr[2]}{$arr[1]}{'ref'}{$arr[7]}++;
}
close IN;


my $client = MongoDB::MongoClient->new(
   host => '10.10.9.35:31112',
   username => 'lintop',
   password => 'lintop.hx321.mongo'
);
$client->reconnect;


my $collection = $client->ns( $db_col_name );


foreach my $hpo(sort {$a cmp $b} keys(%h_hpo2gene)){
	foreach my $gene(sort {$a cmp $b} keys(%{$h_hpo2gene{$hpo}})){
		my %h_record= ();
		$h_record{'gene'} = $gene;
		$h_record{'hpo'} = $hpo;
		$h_record{'match_type'} = $h_hpo2gene{$hpo}{$gene}{'type'};
		my @temp_arr =();
		if(defined($h_hpo2gene{$hpo}{$gene}{'ref'})){
			foreach my $ref(sort {$a cmp $b} keys(%{$h_hpo2gene{$hpo}{$gene}{'ref'}})){
				my %h_temp = ();
				$h_temp{'name'} = $ref;
				push @temp_arr,\%h_temp;
			}
		}
		$h_record{'ref'} = \@temp_arr;
		$collection->insert_one( \%h_record );
	}
}


$client->disconnect;
