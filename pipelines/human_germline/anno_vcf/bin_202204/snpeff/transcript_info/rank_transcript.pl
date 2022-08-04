my %h_hgmd_trans = ();
open HGMD,$ARGV[0] or die $!;
while( $line = <HGMD> ){
	if( $line =~ /Mutation="([^:]+):([^:]+)+:/ ){
		$h_hgmd_trans{$2}++;
	}
	else{
		print STDERR "Illegal hgmd:$line\n";
	}
}
close HGMD;



my %h_all_trans = ();
#my %h_type = ();
my @types = ('HGMD','mRNA','ncRNA','tRNA','rRNA','misc_RNA','precursor_RNA');
open IN,"gunzip -c $ARGV[1]|"or die $!;
while( $line = <IN> ){
	if( $line =~ /\sexon\s(\d+)\s(\d+)\s.*gbkey=([^;]+);gene=([^;]+);.*transcript_id=([^\.;]+)\./){
		my ($pos1,$pos2,$gbkey,$gene,$trans) = ($1,$2,$3,$4,$5);
		if( defined($h_hgmd_trans{$trans}) ){
			$h_all_trans{$gene}{'HGMD'}{$trans} += ($pos2 - $pos1 + 1);
		}
		else{
			$h_all_trans{$gene}{$gbkey}{$trans} += ($pos2 - $pos1 + 1);
		}
	}
	elsif( $line =~ /^NC_012920/ && $line =~ /\sexon\s(\d+)\s(\d+)\s.*gbkey=([^;]+);gene=([^;]+)[;\s]/ ){
		my ($pos1,$pos2,$gbkey,$gene) = ($1,$2,$3,$4);
		foreach my $trans( $gene, "TRANSCRIPT_$gene" ){
			if( defined($h_hgmd_trans{$trans}) ){
				$h_all_trans{$gene}{'HGMD'}{$trans} += ($pos2 - $pos1 + 1);
			}
			else{
				$h_all_trans{$gene}{$gbkey}{$trans} += ($pos2 - $pos1 + 1);
			}
		}
	}
}
close IN;


foreach my $gene( sort {$a cmp $b} keys(%h_all_trans) ){
	print $gene;
	foreach my $type( @types  ){
		foreach my $trans( sort {$h_all_trans{$gene}{$type}{$b} <=> $h_all_trans{$gene}{$type}{$a} } keys(%{$h_all_trans{$gene}{$type}}) ){
			print "\t".$trans."\t".$h_all_trans{$gene}{$type}{$trans};
		}
	}
	print "\n";
}
