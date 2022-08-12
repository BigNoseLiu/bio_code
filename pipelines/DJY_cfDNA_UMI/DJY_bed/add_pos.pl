my %h_pos = ();
open IN,$ARGV[0] or die $!;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_pos{$arr[0]}{'+'} = $arr[1];
	$h_pos{$arr[0]}{'-'} = $arr[2];
}
close IN;

open IN,$ARGV[1] or die $!;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	for(my$i=1;$i<scalar@arr;$i++){
		$arr[$i] =~ s/AGATGTGTATAAGAGACAG//;
		if( defined($h_pos{$arr[$i]}) ){
			$arr[$i] .= "\t".$h_pos{$arr[$i]}{'+'}."\t".$h_pos{$arr[$i]}{'-'};
		}
		else{
			$arr[$i] .= "\t-\t-";
		}
	}
	print join("\t",@arr)."\n";
}
close IN;
