#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use MongoDB;

#my $client = MongoDB->connect( '10.10.9.22:27017' );

my $client = MongoDB::MongoClient->new(
   host => '10.10.9.22:27017'
);
$client->reconnect;

my $collection = $client->ns( $ARGV[0] );


my @data = $collection->find()->all();
foreach my $t(@data){
	my $id =%$t{'_id'};
	$collection->delete_one({'_id'=>$id});
}
$client->disconnect;
