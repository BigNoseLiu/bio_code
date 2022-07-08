my %h_genes = ();
while( $line = <STDIN> ){
	chomp $line;
	$h_genes{$line}{'flag'} = 0;
}

open IN,"gunzip -c $ARGV[0]|" or die $!;
while( $line = <IN> ){
	chomp $line;
	next if($line !~ /^NC_/);
	if( $line =~ /ID=([^;]+);.*Parent=([^;]+);.*;gene=([^;]+);.+$/ ){
		my $exon=$1;
		my $parent=$2;
		my $gene=$3;
		next if($exon !~ /NM_/ && $exon !~ /NP_/ && $gene ne "RMRP" && $gene ne "RNR1");
		my @arr = split(/\t/,$line);
		#if($arr[2] eq "exon" && defined($h_genes{$gene})){
		if(defined($h_genes{$gene})){
			my $chr = $arr[0];
			$chr =~ s/NC_012920.1/chrM/;
			$chr =~ s/^NC_0+([^0]\d*)\.\d*$/chr$1/;
			$chr =~ s/chr23/chrX/;
			$chr =~ s/chr24/chrY/;
			my $pos1 = $arr[3]-1;
			my $pos2 = $arr[4];
			my $out_str = "$chr\t$pos1\t$pos2\t$gene\t$exon;$parent:$chr:$pos1:$pos2";
			if($arr[2] eq "exon"){
				$h_genes{$gene}{'exon'}{$out_str}++;
			}
			if($arr[2] eq "CDS"){
				$h_genes{$gene}{'CDS'}{$out_str}++;
			}
		}
	}
}
close IN;

my $out_dir = $ARGV[1];
`mkdir -p $out_dir/exon $out_dir/CDS`;
foreach my $gene( sort {$a cmp $b} keys(%h_genes) ){
	if( defined($h_genes{$gene}{'exon'})){
		open OUT,">$out_dir/exon/$gene.exon" or die $!;
		foreach my $out_str( keys(%{$h_genes{$gene}{'exon'}}) ){
			print OUT $out_str."\n";
		}
		close OUT;
	}
	else{
		print STDERR "Err: no info for gene exon:$gene\n";
	}

	if( defined($h_genes{$gene}{'CDS'})){
		open OUT,">$out_dir/CDS/$gene.CDS" or die $!;
		foreach my $out_str( keys(%{$h_genes{$gene}{'CDS'}}) ){
			print OUT $out_str."\n";
		}
		close OUT;
	}


}
