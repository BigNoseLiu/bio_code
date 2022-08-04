use LWP::Simple qw/get/;
`mkdir -p dir_citations`;
while( $line = <STDIN> ){
        chomp $line;
        next if( $line !~ /\S/ );
        my $a = get("https://www.pharmgkb.org/clinicalAnnotation/$line");
        if( not -f "dir_citations/pubmed.$line.txt"){
        open OUT,">dir_citations/pubmed.$line.txt" or die $!;
        print OUT $a;
        close OUT;
        }
        sleep(10);
}
