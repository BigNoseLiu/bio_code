#input1: perl $0 sample1.filter.vcf sample2.filter.vcf .......  sampeN.filter.vcf >stat.xls
#input2: perl $0 id1:sample1.filter.vcf id2:sample2.filter.vcf .......  idN:sampeN.filter.vcf >stat.xls
my %h_sex_stat = ();
my @types = ("other","het","hom");
foreach my $data( @ARGV ){

	my ( $file, $file_path ) = ($data,$data);
	if( $data =~ /:/ ){
		($file,$file_path) = split(/:/,$data);
	}
	foreach my $type( @types ){
		$h_sex_stat{$file}{'chrY'}{$type} = 0;
		$h_sex_stat{$file}{'chrX'}{$type} = 0;
	}
	open IN,$file_path or die $!;
	my $line;
	while( $line = <IN> ){
		chomp $line;
		my @arr = split(/\t/,$line);
		next if( !defined($arr[9]) );
		next if( !defined($arr[6]) );
		next if( $arr[6] !~ /PASS/i && $arr[6] ne '.' );
		my $type = $types[0];
		if( $arr[9] =~ /^0\/1:/i || $arr[9] =~ /^1\/0:/i ){
			$type = $types[1];
		}
		elsif( $arr[9] =~ /^1\/1:/i ){
			$type = $types[2];
		}
		if( $arr[0] =~ /Y/i ){
			$h_sex_stat{$file}{'chrY'}{$type}++;
		}
		elsif( $arr[0] =~ /X/i ){
			$h_sex_stat{$file}{'chrX'}{$type}++;
		}
		else{
			next;
		}
	}
	close IN;

}

print join( "\t", 'sample', 'sex_predict', 'chrX', 'chrY');
foreach my $chr('chrX','chrY'){
	foreach my $type( @types ){
		print "\t$chr\_$type";
	}
}
print "\n";
foreach my $file( sort {$a cmp $b} keys( %h_sex_stat )){
	#预测性别
	my $sex = "-";
	if( $h_sex_stat{$file}{"chrX"}{"het"} > 0 && $h_sex_stat{$file}{"chrX"}{"hom"} > 0){
		my $ratio = $h_sex_stat{$file}{"chrX"}{"hom"}/$h_sex_stat{$file}{"chrX"}{"het"};
		if( $ratio > 3 ){
			$sex = "male";
		}
		elsif( $ratio < 1.5 ){
			$sex = "female";
		}
		else{
			$sex = "NA";
		}

	}
	
	my $total_x = $h_sex_stat{$file}{"chrX"}{"hom"}+$h_sex_stat{$file}{"chrX"}{"het"}+$h_sex_stat{$file}{"chrX"}{"other"};
	my $total_y = $h_sex_stat{$file}{"chrY"}{"hom"}+$h_sex_stat{$file}{"chrY"}{"het"}+$h_sex_stat{$file}{"chrY"}{"other"};
	print "$file\t$sex\t$total_x\t$total_y";
	foreach my $chr('chrX','chrY'){
		foreach my $type( @types ){
			print "\t".$h_sex_stat{$file}{$chr}{$type};
		}
	}
	print "\n";
}
