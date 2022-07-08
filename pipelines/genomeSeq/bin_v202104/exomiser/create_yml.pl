while( $line = <STDIN> ){
	foreach my $input( @ARGV ){
		next if( $input !~ /:/ );
		my @arr = split(/:/,$input);
		my $be_replace = shift @arr;
		my $to_replace = join( ":", @arr );
		$line =~ s/$be_replace/$to_replace/;
	}
	print $line;
}
