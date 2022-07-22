while( $line = <STDIN> ){
	$line =~ s/CLNREVSTAT=[^;]*single_submitter[^;]*;/CLNREVSTAT=single_submitter;/;
	$line =~ s/CLNREVSTAT=[^;]*multiple_submitters[^;]*;/CLNREVSTAT=multiple_submitters;/;
	$line =~ s/CLNREVSTAT=[^;]*conflicting_interpretations[^;]*;/CLNREVSTAT=conflicting;/;
	$line =~ s/,/\/\//g;
	print $line;
}
