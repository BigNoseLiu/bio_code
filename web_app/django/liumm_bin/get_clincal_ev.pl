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


open IN,"mut2drug.txt" or die $!;
$line = <IN>;
my %h_mtu2drug2dis = ();
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	next if($arr[0] !~ /\S/);
	$h_mut2drug2dis{$arr[2].$arr[1]} = $arr[0];
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


open OUTEV,">ev.xls" or die $!;
open OUTRAW,">mut2ev.xls" or die $!;
open IN,"clincial_ev.txt" or die $!;
$line = <IN>;
#my %h_drug = ();
my $sum_id = 1000000;
my %h_sum = ();
my %h_evname= ();
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$arr[1] =~ s/rearrange/Rearrange/;
	$arr[1] =~ s/fusion/Fusion/;
	$arr[1] = "ALK Fusion/Rearrange" if( $arr[1] eq "ALK Fusion" );
	$arr[1] = "ROS1 Fusion/Rearrange" if( $arr[1] eq "ROS1 Fusion" );
	my @arr_mut = ($arr[1]);
	if( $arr[3] =~ /\S/ ){
		$arr[2] = $arr[3];
	}
	my @arr_disease = ($arr[2]);
	@arr_disease = ("卵巢癌","输卵管癌","原发性腹膜癌") if($arr[2] eq "卵巢癌包括输卵管癌和原发性腹膜癌");
	@arr_disease = ("食道癌","胃食管交界处癌") if($arr[2] eq "食道癌和胃食管交界处癌");
	my @arr_drug = split(/、/,$arr[9]);
	my $type = $arr[8];
        $type =~ s/敏感/sensitive/;
        $type =~ s/耐药/resistance/;
	my $clin_sum = $arr[10];
	my $clin_name = $arr[11];
	my $ref_level = "";
	if( $clin_sum =~ /FDA/ ){
		$clin_name = "FDA";
		$ref_level = "104_FDA";
	}
	elsif( $clin_sum =~ /NMPA/ ){
		$clin_name = "NMPA";
		$ref_level = "103_NMPA";
	}
	elsif( $clin_sum =~ /NCCN/ ){
		$ref_level = "101_NCCN";
	}
	elsif( $clin_sum =~ /CSCO/ ){
		$ref_level = "102_CSCO";
	}
	else{
		print STDERR $clin_sum."\n";
	}
	if( $clin_name !~ /\S/ ){
		$clin_name = $clin_sum;
	}
	if( !defined($h_sum{$clin_name})){
		$sum_id++;
		$h_sum{$clin_name} = $sum_id;
	}
	$h_evname{$clin_name}{'name'} = "$arr[11]";
	$h_evname{$clin_name}{'version'} = "$arr[12]";


	next if( $type !~ /\S/ );
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
				$t_drug = '塞尔帕替尼' if($t_drug eq 'selpercatinib');
				$t_drug = 'Fam-Trastuzumab Deruxtecan-Nxki' if($t_drug eq 'Fam-trastuzumab deruxtecan-nxki');
				$t_drug = 'Fam-Trastuzumab Deruxtecan-Nxki' if($t_drug eq 'fam-trastuzumab deruxtecan-nxki');
				if( !defined($h_drug{$t_drug}) ){
					print STDERR "$t_drug\n";
				}

				if( !defined($h_drug2dis{$h_drug{$t_drug}.$h_disease{$t_dis}}) ){
					print STDERR $h_drug{$t_drug}."\t".$h_disease{$t_dis}."\n";
					print STDERR $t_drug."\t".$t_dis."\n";
				}
				if( !defined( $h_mut2drug2dis{ $h_mut{$t_mut}.$h_drug2dis{$h_drug{$t_drug}.$h_disease{$t_dis}} }) ){
					print $h_drug2dis{$h_drug{$t_drug}.$h_disease{$t_dis}}."\t".$h_mut{$t_mut}."\n";
					print STDERR $t_drug."\t".$t_dis."\t".$t_mut."\n";
				}
				$arr[6] =~ s/级//;
				print OUTRAW $h_mut2drug2dis{ $h_mut{$t_mut}.$h_drug2dis{$h_drug{$t_drug}.$h_disease{$t_dis}} }."\t$type\tsupport\t$arr[6]\t".$h_sum{$clin_name}."\t$clin_sum\n";
			}
		}
	}

}
close IN;

foreach my $clin_sum( sort {$h_sum{$a} <=> $h_sum{$b}} keys(%h_sum) ){
	print OUTEV $h_sum{$clin_sum}."\t".$clin_sum."\t".$h_evname{$clin_sum}{'name'}."\t".$h_evname{$clin_sum}{'version'}."\n";
}


close OUTEV;
close OUTRAW;
