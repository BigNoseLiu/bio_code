my %h_pos = ();
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	my $id = shift @arr;
	foreach my $pos(split(/,/,join(",",@arr))){
		$h_pos{$pos}{$id}++;
	}
}

foreach my $pos( sort {$a cmp $b} keys(%h_pos) ){
	if( $pos =~ /^(.+):(.+)_(.+)$/ ){
		print $1."\t".($2-100)."\t".($3+100)."\t".join(",",sort {$a cmp $b} keys(%{$h_pos{$pos}}))."\n";
	}
}
