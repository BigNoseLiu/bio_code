my %h_stat = ();
while($file = <STDIN>){
	chomp $file;
	open IN,"gunzip -c $file|" or die $!;
	while($line = <IN>){
		$line = <IN>;
		chomp $line;
		my @arr = split(//,$line);
		foreach my $i(0,1,2,3,4,5,6,7,8,-1,-2,-3,-4,-5,-6,-7,-8,-9){
			$h_stat{$file}{$i}{$arr[$i]}++;
		}
		$line = <IN>;
		$line = <IN>;
	}
	close IN;
}
my @bases = ('A','T','C','G','N');
print "#file\trow\t".join("\t",@bases)."\n";
foreach my $file(sort {$a cmp $b} keys(%h_stat)){
	foreach my $row(sort {$a cmp $b} keys(%{$h_stat{$file}})){
		print "$file\t$row";
		foreach my $base( @bases ){
			if( defined( $h_stat{$file}{$row}{$base} ) ){
				print "\t".$h_stat{$file}{$row}{$base};
			}
			else{
				print "\t0";
			}
		}
		print "\n";
	}
}
