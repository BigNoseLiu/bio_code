my %h_pos = ();
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	my $id = shift @arr;
	$id =~ s/.\d$//;
	$id =~ s/-$//;
	foreach my $pos(split(/,/,join(",",@arr))){
		if( $pos =~ /^(.+):(.+)_(.+)$/ ){
			my ($chr,$pos1,$pos2) = ($1,$2,$3);
			$h_pos{$id}{'chr'} = $chr;
			$h_pos{$id}{'pos1'} = $pos1 if(!defined($h_pos{$id}{'pos1'})|| $h_pos{$id}{'pos1'} > $pos1);
			$h_pos{$id}{'pos2'} = $pos2 if(!defined($h_pos{$id}{'pos2'})|| $h_pos{$id}{'pos2'} < $pos2);
		}
	}
}

foreach my $id( sort {$a cmp $b} keys(%h_pos) ){
	print $h_pos{$id}{'chr'}."\t".$h_pos{$id}{'pos1'}."\t".$h_pos{$id}{'pos2'}."\t".$id."\n";
}
