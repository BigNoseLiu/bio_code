while( $line = <STDIN> ){
	if( $line =~ /\dinv/ ){
		print $line;
	}
}
