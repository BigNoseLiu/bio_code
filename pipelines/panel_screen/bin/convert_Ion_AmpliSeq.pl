$line = <STDIN>;
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	if($arr[7] =~ /GENE_ID=([^;]+);/){
		my $gene = $1;
		print "$arr[0]\t$arr[1]\t$arr[2]\t$gene\t$arr[7];ID2=$arr[3]\n";
	}
}
