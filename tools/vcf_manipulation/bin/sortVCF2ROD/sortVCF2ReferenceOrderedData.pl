use strict;
open FAI,$ARGV[0] or die $!;
my @lines = <FAI>;
#my %h_order = ();
my $count = 0;
my @orders = ();
foreach my $line ( @lines ){
	chomp $line;
	my @arr = split(/\s+/,$line);
	#$count++;
	#$h_order{ $arr[0] } = $count;
	push @orders,$arr[0];
}
close FAI;

@lines = <STDIN>;
my $contig_pos = 0;
my %h_contigs = ();
my %h_variants = ();
foreach my $line ( @lines ){
	if( $line =~ /^##contig.*ID=([^,]+),/){
		$h_contigs{$1} = $line;
		$contig_pos = 1;
		next;
	}

	if( $line=~ /^#/ && $line !~ /^##contig/){
		if( $contig_pos == 1 ){
			foreach my$chr( @orders ){
				if( defined( $h_contigs{$chr} ) ){
					print $h_contigs{$chr};
				}
				else{
					print STDERR "Err: unmatch ref & variant chrs\t$chr\n";
					exit(0);
				}
			}
		}

		$contig_pos = 0;
		print $line;
		next;
	}


	chomp $line;
	my @arr = split(/\s+/,$line);
	$h_variants{ $arr[0] }{$arr[1]} = $line;
}


foreach my$chr( @orders ){
	if( defined( $h_variants{$chr} ) ){
		foreach my $pos( sort {$a <=> $b } keys(%{$h_variants{$chr}}) ){
			print $h_variants{$chr}{ $pos }."\n";
		}
	}
}
