print "#Chr\tStart\tEnd\tRef\tAlt\tmitofreq.GnomAD_ALL.freq\tmitofreq.GnomAD_ALL.AF_hom\tmitofreq.GnomAD_ALL.hom\tmitofreq.GnomAD_ALL.AF_het\tmitofreq.GnomAD_ALL.het\tmitofreq.GnomAD_ALL.max_het\n";
$line = <STDIN>;
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	my ($pos,$ref,$alt) = ($arr[1],$arr[2],$arr[3]);
	my ($pos1,$pos2) = ($pos,$pos+length($ref)-1);
	if( $alt =~ s/^$ref(\S+)$/$1/ ){
		$pos1 = $pos + length($ref) - 1 ;
		$pos2 = $pos + length($ref) - 1 ;
		$ref = "-";
	}
	elsif( $ref =~ s/^(\S+)$alt$/$1/ ){
		$pos1 = $pos;
		$pos2 = $pos + length($ref) - 1;
		$alt = "-";
	}
	elsif( $ref =~ s/^$alt(\S+)$/$1/ ){
		$pos1 = $pos + length($alt);
		$alt = "-";
	}
	for my $i(7,8,10){
		$arr[$i] = sprintf("%f",$arr[$i]);
		$arr[$i] =~ s/([^0])0*$/$1/;
		$arr[$i] =~ s/^0\.$/0/;
	}
	print "MT\t$pos1\t$pos2\t$ref\t$alt\t".($arr[7]+$arr[8])."\t$arr[7]\t$arr[5]\t$arr[8]\t$arr[6]\t$arr[10]\n";
}
