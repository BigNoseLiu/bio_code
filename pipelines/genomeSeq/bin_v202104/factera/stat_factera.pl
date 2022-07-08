#input1: perl $0 sample1.factera.fusions.txt sample2.factera.fusions.txt ..... sampleN.factera.fusions.txt >factera.fusions.stat.xls
#input2: perl $0 id1:sample1.factera.fusions.txt id2:sample2.factera.fusions.txt ..... idN:sampleN.factera.fusions.txt >factera.fusions.stat.xls
sub mysort{

	my @arr_a = split(/\t/,$a);
	my @arr_b = split(/\t/,$b);
	if( $arr_a[1] gt $arr_a[2] ){
		my $t = $arr_a[1];
		$arr_a[1] = $arr_a[2];
		$arr_a[2] = $t;
	}

	if( $arr_b[1] gt $arr_b[2] ){
		my $t = $arr_b[1];
		$arr_b[1] = $arr_b[2];
		$arr_b[2] = $t;
	}
	return join("\t",@arr_a) cmp join("\t",@arr_b);
}

my $header = "";
my %h_factera = ();
my @samples = ();
foreach my $data( @ARGV ){

	my ( $file, $file_path ) = ($data,$data);
	if( $data =~ /:/ ){
		($file,$file_path) = split(/:/,$data);
	}
	push @samples,$file;


	open IN,$file_path or die $!;
	my $line;
	$header = <IN>;
	while( $line = <IN> ){
		chomp $line;
		my @arr = split(/\t/,$line);
		my $tag = "PASS";
		my $ratio1 = $arr[6]/$arr[5];
		my $ratio2 = $arr[5]/$arr[6];
		if( $arr[5] < 5 || $arr[6] < 5 ){
			$tag = "FILTER_LowDepth";
		}
		elsif( $ratio1 > 3 || $ratio2 > 3 ){
			$tag = "FILTER_HighRatio";
		}
		$h_factera{join("\t",$arr[0],$arr[1],$arr[2])}{$file}{join(":",$tag,$arr[3],$arr[4],$arr[5],$arr[6])}{"tag"} = $tag; 
	}
	close IN;

}
print join("\t","Est_Type\tRegion1\tRegion2\tFilter\tPair",@samples)."\n";
#foreach my $fusion( sort {$a cmp $b} keys(%h_factera)){
foreach my $fusion( sort mysort keys(%h_factera)){
	#check pair status
	my $pair = "No_Pair";
	my @arr= split(/\t/,$fusion);
	my $rev_fusion = $arr[0]."\t".$arr[2]."\t".$arr[1];
	foreach my $file( keys(%{$h_factera{$fusion}})){
		if( defined( $h_factera{$rev_fusion}{$file} ) ){
			$pair = "Pair";
		}
	}
	my $tag = "FILTER";
	my $detail_str = "";
	foreach my $file( @samples ){
		if( defined( $h_factera{$fusion}{$file} ) ){
			my @infos = sort {$a cmp $b} keys(%{$h_factera{$fusion}{$file}});
			foreach my $info( @infos ){
				if( $h_factera{$fusion}{$file}{$info}{"tag"} eq "PASS" ){
					$tag = "PASS";
				}
			}
			$detail_str .= "\t".join("&", @infos );
		}
		else{
			$detail_str .= "\t-";
		}
	}
	print "$fusion\t$tag\t$pair$detail_str\n";
}

