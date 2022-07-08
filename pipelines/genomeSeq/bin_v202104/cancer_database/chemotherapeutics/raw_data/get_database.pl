my %h_rs_info = ();
open IN,$ARGV[0] or die $!;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	#3	14187449	14187449	G	T	rs2228001
	$h_rs_info{lc($arr[5])}{join(":",@arr[0..4])}++;
}
close IN;


print "分类	药物(中文)	药物（英文）	基因	证据等级	rs号	变异信息	变异结果	解读结果\n";
open IN,$ARGV[1] or die $!;
my $line = <IN>;
my $class = "-";
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	if( $arr[0] =~ /\S/ ){
		$class = $arr[0];
	}
	my $drug_ch = $arr[1];
	my $drug_en = $arr[1];
	if( $arr[1] =~ /^\s*(.+)\((.+)\)\s*$/ ){
		$drug_en = $1;
		$drug_ch = $2;
	}
	my $rs_info = "-";
	if( !defined($h_rs_info{lc($arr[3])}) ){
		print STDERR "$arr[3]\n";
	}
	my $rs_info = lc($arr[3])."\t".join('|',sort {$a cmp $b} keys(%{$h_rs_info{lc($arr[3])}}));
	print "$class\t$drug_ch\t$drug_en\t$arr[2]\t$arr[10]\t$rs_info\t$arr[4]\t$arr[5]\n";
	print "$class\t$drug_ch\t$drug_en\t$arr[2]\t$arr[10]\t$rs_info\t$arr[6]\t$arr[7]\n";
	print "$class\t$drug_ch\t$drug_en\t$arr[2]\t$arr[10]\t$rs_info\t$arr[8]\t$arr[9]\n";
}
close IN;
