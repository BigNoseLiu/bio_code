#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use MongoDB;

#my $client = MongoDB->connect( '10.10.9.22:27017' );

#   host => '1.12.236.215:27017',
#   host => 'gmbzero.tpddns.cn:31112',
my $client = MongoDB::MongoClient->new(
   host => '10.10.9.35:31112',
   username => 'lintop',
   password => 'lintop.hx321.mongo'
);
$client->connect;



my $collection = $client->ns( $ARGV[0] );


my @data = $collection->find()->all();
foreach my $t(@data){
	my $id =%$t{'_id'};
	$collection->delete_one({'_id'=>$id});
}
$client->disconnect;
