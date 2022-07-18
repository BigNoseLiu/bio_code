while( $line = <STDIN> ){
	$line =~ s/,/\/\//g;
	print $line;
}
