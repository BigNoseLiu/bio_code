#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use MongoDB;

my $client = MongoDB::MongoClient->new(
   host => '1.12.236.215:27017',
   username => 'lintop',
   password => 'lintop.hx321.mongo'
);
$client->reconnect;




my $collection = $client->ns( $ARGV[0] );


my @data = $collection->find()->all();
foreach my $t(@data){
	my $id =%$t{'_id'};
	$collection->delete_one({'_id'=>$id});
}
$client->disconnect;
