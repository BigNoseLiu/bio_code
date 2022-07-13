my @arr_hand = ();
my %h_header = ();
my %h_lines = ();
my $f_count = 0;
while( $file = <STDIN> ){
	chomp $file;
	open $arr_hand[$f_count],$file or die $!;
	$t = $arr_hand[$f_count];
	$header = <$t>;
	chomp $header;
	$h_header{$f_count} = $header;
	$f_count++;
}

foreach my $f_hand(@arr_hand){
	close $f_hand;
}
