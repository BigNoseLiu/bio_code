$line = <STDIN>;
print "#Chr\tStart\tEnd\tRef\tAlt\tMitoTIP_score\n";
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$arr[2] =~ s/:/-/;
	print "MT\t$arr[0]\t$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\n";
}
