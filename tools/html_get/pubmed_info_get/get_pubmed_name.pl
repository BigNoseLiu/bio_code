use strict;
use LWP::UserAgent;
my $get_page = LWP::UserAgent -> new;
my $html = $get_page -> get("https://www.ncbi.nlm.nih.gov/pubmed/17606709");
my $line;
while( $line = <STDIN> ){
	next if( $line =~ // );
}
print $html->content."\n";
