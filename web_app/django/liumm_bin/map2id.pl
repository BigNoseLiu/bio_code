open IN,$ARGV[0] or die $!;
my $col1 = $ARGV[1];
$line = <IN>;
my %h_id = ();
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	next if($arr[0] !~ /\S/);
	$h_id{$arr[$col1]} = $arr[0];
}
close IN;



open IN,$ARGV[2] or die $!;
while( $line = <IN> ){
	chomp $line;
	if(defined($h_id{$line})){
		print $h_id{$line}."\n";
	}
	else{
		print "$line\n";
	}
}
close IN;
