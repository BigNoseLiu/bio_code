my $header = "";

my %h_cnv = ();
foreach my $file( @ARGV ){

	open IN,$file or die $!;
	my $line;
	$header = <IN>;
	while( $line = <IN> ){
		chomp $line;
		my @arr = split(/\t/,$line);
		$h_cnv{join("\t",$arr[0],$arr[1],$arr[2],$arr[3])}{$file}{join(":",$arr[4])}++; 
	}
	close IN;

}

print join("\t","chromosome\tstart\tend\tgene",@ARGV)."\n";
foreach my $cnv( sort {$a cmp $b} keys(%h_cnv)){
	print "$cnv";
	foreach my $file( @ARGV ){
		if( defined($h_cnv{$cnv}{$file}) ){
			print "\t".join("&", sort {$a cmp $b} keys(%{$h_cnv{$cnv}{$file}}) );
		}
		else{
			print "\t-";
		}
	}
	print "\n";
}
