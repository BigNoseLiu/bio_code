while( $line = <STDIN> ){
	next if($line =~ /^>/);
	next if($line !~ /\S/);
	chomp $line;
	my @arr = split(/\t/,$line);
	my $rs = "rs$arr[4]";
	my $chr = "chr$arr[0]";
	my $pos1 = $arr[1];
	my $gene = $arr[8];
	if( $arr[2] =~ /^([ATCG-]+)>(.*)$/ ){
		my $ref = $1;
		my $pos2 = $pos1 + length($ref) -1;
		foreach my $alt( split(/,/,$2) ){
			print "$rs\t$gene\t$rs\t$chr\t$pos1\t$pos2\t$ref\t$alt\t$arr[3]\n";
		}
	}
	else{
		print STDERR $line."\n";
	}
}
#CYP2C19*1       CYP2C19 rs3758581       NC_000010.11    94842866        94842866        A       G       substitution
