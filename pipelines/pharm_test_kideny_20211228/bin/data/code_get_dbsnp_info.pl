my $query_str = "";
my $count = 0;
while( $line = <STDIN> ){
	$count++;
	chomp $line;
	my @arr = split(/\t/,$line);
	$query_str .= "or $arr[1]".'[Reference SNP ID]';
}
$query_str =~ s/^or//;
my $cmd = "esearch -db snp -query '".$query_str."' |efetch -format tab";
my $cmd_str = `$cmd`;
print ">$cmd\n";
print $cmd_str;
