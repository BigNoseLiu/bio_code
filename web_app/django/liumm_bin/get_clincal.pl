open IN,"mut.txt" or die $!;
$line = <IN>;
my %h_mut = ();
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	next if($arr[0] !~ /\S/);
	$h_mut{$arr[3]} = $arr[0];
}
close IN;


open IN,"disease.txt" or die $!;
$line = <IN>;
my %h_disease = ();
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	next if($arr[0] !~ /\S/);
	$h_disease{$arr[3]} = $arr[0];
}
close IN;


open IN,"drug.txt" or die $!;
$line = <IN>;
my %h_drug = ();
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	next if($arr[0] !~ /\S/);
	$h_drug{$arr[3]} = $arr[0];
}
close IN;




open IN,"drug2disease.txt" or die $!;
$line = <IN>;
my %h_drug2dis = ();
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	next if($arr[0] !~ /\S/);
	$h_drug2dis{$arr[2].$arr[1]} = $arr[0];
}
close IN;


open IN,"clincial_info.txt" or die $!;
$line = <IN>;
#my %h_drug = ();
my $sum_id = 100000;
my %h_sum = ();
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$arr[1] =~ s/rearrange/Rearrange/;
	my @arr_mut = ($arr[1]);
	my @arr_disease = ($arr[2]);
	@arr_disease = ("卵巢癌","输卵管癌","原发性腹膜癌") if($arr[2] eq "卵巢癌包括输卵管癌和原发性腹膜癌");
	@arr_disease = ("食道癌","胃食管交界处癌") if($arr[2] eq "食道癌和胃食管交界处癌");
	my @arr_drug = split(/、/,$arr[6]);
	my $type = $arr[7];
        $type =~ s/敏感/sensitive/;
        $type =~ s/耐药/resistance/;
	my $clin_sum = $arr[8];
	if(!defined($h_sum{$clin_sum."\t".$type})){
		$sum_id++;
		$h_sum{$clin_sum."\t".$type} = $sum_id;
	}

	foreach my $t_mut( @arr_mut ){
		if( !defined($h_mut{$t_mut}) ){
			print STDERR "$t_mut\n";
		}
		foreach my $t_dis( @arr_disease){
			if( !defined($h_disease{$t_dis}) ){
				print STDERR "$t_dis\n";
			}
			foreach my $t_drug( @arr_drug){
				$t_drug = '瑞派替尼' if($t_drug eq 'ripretinib');
				$t_drug = '舒尼替尼' if($t_drug eq 'sunitinib');
				$t_drug = '瑞戈非尼' if($t_drug eq 'regorafenib');
				$t_drug = '塞尔帕替尼' if($t_drug eq 'Selpercatinib');
				$t_drug = '艾伏尼布' if($t_drug eq 'ivosidenib');
				$t_drug = '帕尼单抗' if($t_drug eq 'Panitumumab');
				$t_drug = '特泊替尼' if($t_drug eq 'Tepotinib');
				$t_drug = '卡马替尼' if($t_drug eq 'Capmatinib');
				$t_drug = '恩曲替尼' if($t_drug eq 'entrectinib');
				$t_drug = '芦卡帕利' if($t_drug eq 'rucaparib');
				$t_drug = '莫博赛替尼' if($t_drug eq 'Mobocertinib');
				$t_drug = '英菲格拉替尼' if($t_drug eq 'infigratinib');
				if( !defined($h_drug{$t_drug}) ){
					print STDERR "$t_drug\n";
				}
				if( !defined($h_drug2dis{$h_drug{$t_drug}.$h_disease{$t_dis}}) ){
					print STDERR $h_drug{$t_drug}."\t".$h_disease{$t_dis}."\n";
				}
				print $h_mut{$t_mut}."\t".$h_drug2dis{$h_drug{$t_drug}.$h_disease{$t_dis}}."\t".$h_sum{$clin_sum."\t".$type}."\n";
			}
		}
	}

}
close IN;

foreach my $clin_sum( sort {$a cmp $b} keys(%h_sum) ){
	#print $h_sum{$clin_sum}."\t".$clin_sum."\n";
}


