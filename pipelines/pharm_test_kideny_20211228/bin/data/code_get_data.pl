my $pharvar_dir = "/home/liumm/biopipe/biodata/database/pharmvar/pharmvar-5.1.5";
my $pharm_dir = "/home/liumm/biopipe/biodata/database/pharmgkb";
my $line = "";

my %h_change_var_name = ();
open IN,"data_changeVarName.txt" or die $!;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_change_var_name{$arr[0]} = $arr[1];
}
close IN;

my %h_filter_var= ();
open IN,"data_filter_variant.txt" or die $!;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_filter_var{$arr[0]}++;
}
close IN;


my %h_mut_freq = ();
my $freq_header = "";
if( -f "data_mut.hg38_multianno.txt" ){
	open IN,"data_mut.hg38_multianno.txt" or die $!;
	$freq_header = <IN>;
	while( $line = <IN> ){
		chomp $line;
		my @arr = split(/\t/,$line);
		$arr[13] = '-1' if( $arr[13] eq '.' );
		$h_mut_freq{ $arr[0]."\t".$arr[1]."\t".$arr[2]."\t".$arr[3]."\t".$arr[4] }{"freq"} = $arr[13];
		$h_mut_freq{ $arr[0]."\t".$arr[1]."\t".$arr[2]."\t".$arr[3]."\t".$arr[4] }{"all"}{$line}++;
	}
	close IN;
}



my %h_hap = ();
open OUT,">$ARGV[1].freq.xls" or die $!;
print OUT "#var\t$freq_header";
open IN,"cat data_rs.GRCh38.haplotypes.tsv $pharvar_dir/*/GRCh38/*.haplotypes.tsv $pharm_dir/allele_definition_table_202109/GRCh38/*.haplotypes.tsv|" or die $!;
while( $line = <IN> ){
	chomp $line;
	$line =~ s/NC_000023.11/chrX/;
	$line =~ s/NC_000024.10/chrY/;
	$line =~ s/NC_0+([^0]\d*).\d+(\s)/chr$1$2/;
	my @arr = split(/\t/,$line);
	if(scalar@arr > 7 && $arr[8] eq "insertion" && $arr[6] eq "-"){
		$arr[5] = $arr[4];
		$line = join("\t",@arr);
	}
	my $var = $arr[0];
	$var = $h_change_var_name{$var} if( defined($h_change_var_name{$var}) );
	$h_hap{$var}{"def"}{$line}++;
	if( defined($h_mut_freq{ $arr[3]."\t".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[7] }) ){
		$h_hap{$var}{"freq"}{ $h_mut_freq{ $arr[3]."\t".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[7] }{"freq"} }++;
		foreach my $t_freq(keys(%{ $h_mut_freq{ $arr[3]."\t".$arr[4]."\t".$arr[5]."\t".$arr[6]."\t".$arr[7] }{"all"} })){
			print OUT $var."\t".$t_freq."\n";
		}
	}
	else{
		$h_hap{$var}{"freq"}{ '-2' }++;
	}
	#得到ref序列信息
	if( $var =~ /\*/ ){
		$h_hap{$var}{"ref"}{'*1'}++;
	}
	else{
		$h_hap{$var}{"ref"}{$arr[6]}++;
	}
}
close IN;
close OUT;


my %h_target_drug = ();
open IN,$ARGV[0] or die $!;
$line = <IN>;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_target_drug{lc($arr[2])}{"info"} = join("|",@arr);
	$h_target_drug{lc($arr[2])}{"flag"} = 0;
}
close IN;

my %h_ann_alleles = ();
open IN,"$pharm_dir/download_20211228/clinical_ann_alleles.tsv" or die $!;
$line = <IN>;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_ann_alleles{$arr[0]}{ join("\t",$arr[1],$arr[2],$arr[3],$arr[4]) }++;
}
close IN;

my %h_ann_ev= ();
open IN,"$pharm_dir/download_20211228/clinical_ann_evidence.tsv" or die $!;
$line = <IN>;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$arr[-1] = -100 if($arr[-1] =~ /This annotation is not used for clinical annotation scoring/);
	$h_ann_ev{$arr[0]}{ $arr[-1] }++;
}
close IN;


my %h_all_variants = (
'rs15561'=>1,
'rs1057126'=>1,
'rs58973490'=>1,
'rs3758581'=>1,
'rs17879685'=>1
);
open IN,"$pharm_dir/download_20211228/clinical_annotations.tsv" or die $!;
my %h_header = ();
$line = <IN>;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	my $all_drug = $arr[10];
	$all_drug = "compound sulfamethoxazole tablets(sulfamethoxazole/trimethoprim)" if($all_drug =~ /sulfamethoxazole/);
	$all_drug = "felodipine" if($all_drug =~ /calcium channel blockers/);
	foreach my $drug( split(/;/,$all_drug) ){
		if( defined($h_target_drug{lc($drug)}) ){
			#处理变异
			$arr[1] =~ s/Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham/liumm_special_G6PD1/g;
			$arr[1] =~ s/Canton, Taiwan-Hakka, Gifu-like, Agrigento-like/liumm_special_G6PD2/g;
			my $var_list = "";
			my %h_temp_filter_var = ();
			foreach my $var( split(/,/,$arr[1]) ){
				$var =~ s/^\s*(\S)/$1/;
				$var =~ s/(\S)\s*$/$1/;
				if( defined($h_filter_var{$var}) ){
					my $var_simple = $var;
					$var_simple =~ s/^[^\*]+(\*)/$1/;
					$h_temp_filter_var{$var_simple}++;
					next;
				}
				$var =~ s/liumm_special_G6PD1/Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham/g;
				$var =~ s/liumm_special_G6PD2/Canton, Taiwan-Hakka, Gifu-like, Agrigento-like/g;
				$var_list .=";;$var";
			}
			$var_list =~ s/^;;//;
			next if( $var_list =~ /^\s*$/ );
			foreach my $var( split(/;;/,$var_list) ){
				$h_all_variants{$var}++;
			}
			foreach $anno_allele(keys(%{$h_ann_alleles{$arr[0]}})){
				my @arr_annos = split("\t",$anno_allele);
				my $filter_flag = 0;
				foreach my $var_filter( keys(%h_temp_filter_var) ){
					$filter_flag = -1 if($filter_flag!=1 && $filter_flag!=2); 

					my $var_filter_match  = $var_filter;
					$var_filter_match  =~ s/^\*/\\*/;
					#包含无法分型变异的临床说明
					$filter_flag = 1 if($arr_annos[0]." " =~ /$var_filter_match\D/ && $filter_flag != 2);
					#无法分型变异的临床说明
					$filter_flag = 2 if($arr_annos[0] eq $var_filter);
				}
				#去除无法分型的变异
				next if($filter_flag >0);
				#my @support_scores = sort {$b<=>$a} keys(%{$h_ann_ev{$arr[0]}});
				$h_target_drug{lc($drug)}{"var"}{join("\t",$arr[3],$arr[0],$arr[1],$arr[2],$arr[6],$arr[7],$arr[10])}{$anno_allele}++;
				#print STDERR $arr[3]."\n";
				#$h_target_drug{lc($drug)}{"flag"} = 1;
			}
		}
	}
}
close IN;

open OUT,">$ARGV[1].detail.xls" or die $!;
print OUT "#药物\t证据等级\t编号\t变异\t基因\t证据分值\t影响\t药物\tmax频率\tall频率\t类型\t代谢\t疗效\t毒副\t药量\t基因型\t描述\n";
foreach my $drug( sort {$a cmp $b} keys(%h_target_drug)){
	
	#if( $h_target_drug{$drug}{"flag"} == 0){
	#	print STDERR "Err: no info for drug $drug\n";
	#	next;
	#}

	my $str_drug = $h_target_drug{$drug}{"info"};
	my $lev_1_2 = 0;
	my $lev_3 = 0;
	foreach my $var( sort {$a cmp $b} keys(%{$h_target_drug{$drug}{"var"}})){
		my @arr_var = split(/\t/,$var);

		my $ref_base = "NA";
		if( defined( $h_hap{$arr_var[2]}{"ref"} ) ){
			my @arr_refs = keys(%{$h_hap{$arr_var[2]}{"ref"}});
			if( @arr_refs >1 ){
				print STDERR "Waring: multiple ref for $arr_var[2]\n";
			}
			else{
				$ref_base = $arr_refs[0];
			}
		}
		#获取变异频率
		my $max_freq = "-";
		my $min_freq = "-";
		if( defined( $h_hap{$arr_var[2]}{"freq"} ) ){
			my @arr_freqs = sort {$b <=> $a} keys(%{$h_hap{$arr_var[2]}{"freq"}});
			$max_freq = $arr_freqs[0];
			$min_freq = join(",",@arr_freqs);
		}
		next if( $max_freq ne "-" && ($max_freq <0.01||$max_freq >0.99) );
		my $str_var = "$var\t$max_freq\t$min_freq";

		foreach my $allele( sort {$a cmp $b} keys(%{$h_target_drug{$drug}{"var"}{$var}})){

			my @arr_allele = split(/\t/,$allele);

			#只保留高等级的变异
			$lev_1_2 = 1 if($var =~ /^[12]/);
			if($var =~ /^3/){
				$lev_3 = 1;
				next if( $lev_1_2 == 1 );
			}
			if($var =~ /^4/){
				next if( $lev_1_2 == 1 || $lev_3 == 1 );
			}

			#
			my $ref_type = "-";
			if( $arr_allele[0] eq "$ref_base$ref_base" ){
				$ref_type = "wild_type";
			}

			#自动评估疗效等
			next if( $arr_var[4] eq 'Other' );
			my ($toxi,$dosage,$eff,$meta) = ("-","-","-","-");
			if( $arr_var[4] =~ /Toxicity/ ){
				$toxi = "NA";
				if( $allele =~ /increased risk/ || $allele =~ /increased, but not absent, risk/ ){
					$toxi = "high"
				}
				elsif( $allele =~ /decreased risk/ || $allele =~ /decreased, but not absent, risk/ ){
					$toxi = "low"
				}
			}

			if( $arr_var[4] =~ /Efficacy/ ){
				$eff = "NA";
				if( $allele =~ /increased response/ || $allele =~ /better response/ ){
					$eff = "high"
				}
				elsif( $allele =~ /decreased response/ || $allele =~ /worse response/ || $allele =~ /poorer response/ ){
					$eff = "low"
				}
			}
			if( $arr_var[4] eq 'Dosage' ){
				$dosage = "NA";
			}


			if( $arr_var[4] eq 'Metabolism/PK' ){
				$meta= "NA";
				if( $allele =~ /increased/ ){
					$meta = "high"
				}
				elsif( $allele =~ /decreased/ ){
					$meta = "low"
				}
			}
			elsif( $arr_var[4] =~ /Metabolism/ ){
				$meta= "NA";
				if( $allele =~ /increased metabolism/ ){
					$meta = "high"
				}
				elsif( $allele =~ /decreased metabolism/ ){
					$meta = "low"
				}
			}


			print OUT "$str_drug\t$str_var\t$ref_type\t$meta\t$eff\t$toxi\t$dosage\t$allele\n";
			$h_target_drug{lc($drug)}{"flag"} = 1;
			$str_var = "\t\t\t\t\t\t\t\t";
			$str_drug = "";
		}
	}
}
close OUT;


foreach my $drug( sort {$a cmp $b} keys(%h_target_drug)){
	
	if( $h_target_drug{$drug}{"flag"} == 0){
		print STDERR "Err: no info for drug $drug\n";
	}
}

open OUT,">$ARGV[1].var.xls" or die $!;
foreach my $var( sort {$a cmp $b} keys(%h_all_variants)){

	if(defined($h_hap{$var})){
		foreach my $mut( sort {$a cmp $b} keys(%{$h_hap{$var}{"def"}})){
			print OUT $var."\t".$mut."\n";
		}
	}
	else{
		print STDERR $var."\t0\n";
	}
}
close OUT;

#clinical_ann_alleles.tsv
#clinical_ann_evidence.tsv
#clinical_annotations.tsv
#var_drug_ann.tsv

