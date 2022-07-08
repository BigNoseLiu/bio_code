my %h_bed = ();
while( $line = <STDIN> ){
	chomp $line;
	if( $line =~ /^#/ ){
		print $line."\n";
	}
	my @arr = split(/\t/,$line);
	$arr[0] =~ s/^chr//i;
	$arr[0] = "MT" if($arr[0] =~ /m/i);
	my $chr = $arr[0];
	my $pos1 = $arr[1];
	my $pos2 = $arr[2];
	$chr = 23 if( $chr =~ /x/i );
	$chr = 24 if( $chr =~ /y/i );
	$chr = 25 if( $chr =~ /m/i );
	$h_bed{$chr}{$pos1}{$pos2} = join("\t",@arr);
	my $count = scalar@arr;
	if( $count < 4  ){
		$h_bed{$chr}{$pos1}{$pos2} .= "\t$chr-$pos1-$pos2";
	}
	if( $count < 5  ){
		$h_bed{$chr}{$pos1}{$pos2} .= "\t1";
	}
	if( $count < 6  ){
		$h_bed{$chr}{$pos1}{$pos2} .= "\t.";
	}
}

foreach my $chr( sort {$a <=> $b} keys(%h_bed) ){
	foreach my $pos1( sort {$a <=> $b} keys(%{$h_bed{$chr}}) ){
		foreach my $pos2( sort {$a <=> $b} keys(%{$h_bed{$chr}{$pos1}}) ){
			print $h_bed{$chr}{$pos1}{$pos2}."\n";
		}
	}
}
