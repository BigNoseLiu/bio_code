while( $line = <STDIN> ){
	@arr = split(/\t/,$line);
	next if( $arr[2] =~ /_/ );
	$arr[2] =~ s/chr//;
	$arr[2] =~ s/M/MT/;
	print join("\t",@arr);
}

