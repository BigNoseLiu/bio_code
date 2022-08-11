my @out_dirs = `ls out_data`;
chomp @out_dirs;
my %h_out = ();
foreach my $out( @out_dirs ){
	$h_out{$out}++;
}
while($file=<STDIN>){
	chomp $file;
	open IN,$file or die $!;
	print $file."\n";
	while( $line = <IN> ){
		my @arr = split(/\s+/,$line);
		if( !defined($h_out{$arr[4]}) ){
			print $line;
		}
		else{
			if( -f 'out_data/'.$arr[4].'/'.$arr[4].'.09.somatic_raw.vcf'){
			}
			else{
				print STDERR "no raw vcf for $arr[4]\n";
			}
			if( -f 'out_data/'.$arr[4].'/'.$arr[4].'.10.somatic_raw.annovar.hg19_multianno.vcf'){
			}
			else{
				print STDERR "no anno vcf for $arr[4]\n";
			}
		}
	}
	close IN;
}
