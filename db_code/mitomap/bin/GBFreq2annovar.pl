print "#Chr\tStart\tEnd\tRef\tAlt\tmitofreq.GBFreq_ALL.freq\tmitofreq.GBFreq_ALL.AC_all\tmitofreq.GBFreq_African.freq\tmitofreq.GBFreq_Asian.freq\tmitofreq.GBFreq_Eurasian.freq\n";
while( $line = <STDIN> ){
	if( $line !~ /^\d/ ){
		next;
	}
	chomp $line;
	my @arr = split(/\t/,$line);
	$arr[3] =~ s/d/-/;
	my $ref = $arr[1];
	if( $arr[3] =~ s/^$ref(\S+)$/$1/ ){
		$arr[2] = $arr[2] + length($arr[1]) - 1 ;
		$arr[1] = "-";
	}
	for my $i(5,6,7,8){
		$arr[$i] =~ s/%//;
		$arr[$i] = $arr[$i]/100;
	}
	print "MT\t$arr[2]\t$arr[2]\t$arr[1]\t$arr[3]\t$arr[5]\t$arr[4]\t$arr[6]\t$arr[7]\t$arr[8]\n";
	if( $arr[1] =~ /[^ATCG]/ ){
		#print STDERR $arr[1]."\n";
	}
	if( $arr[3] =~ /[^ATCG]/ ){
		#print STDERR $arr[3]."\n";
	}
}
