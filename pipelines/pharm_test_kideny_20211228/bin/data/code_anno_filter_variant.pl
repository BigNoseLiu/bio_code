my $line = "";

my %h_freq = ();
open IN,$ARGV[0] or die $!;
while( $line = <IN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_freq{$arr[0]}{$arr[1]}{$arr[2]}{$arr[3]}{$arr[4]} = $line;
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


my %h_hap = ();
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
	$h_hap{$var}{$line}++;
}
close IN;



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
			foreach $anno_allele(keys(%{$h_ann_alleles{$arr[0]}})){
				my @arr_annos = split("\t",$anno_allele);
				my $filter_flag = 0;
				foreach my $var_filter( keys(%h_temp_filter_var) ){
					$filter_flag = -1 if($filter_flag!=1 && $filter_flag!=2); 

					my $var_filter_match  = $var_filter;
					$var_filter_match  =~ s/^\*/\\*/;
					#包含无法分型变异的临床说明
					$filter_flag = 1 if($arr_annos[0] =~ /$var_filter_match/ && $filter_flag != 2);
					#无法分型变异的临床说明
					$filter_flag = 2 if($arr_annos[0] eq $var_filter);
				}
				#去除无法分型的变异
				next if($filter_flag >0);
				$h_target_drug{lc($drug)}{"var"}{join("\t",$arr[3],$arr[0],$arr[1],$arr[2],$arr[7],$arr[10])}{$anno_allele}++;
				#print STDERR $arr[3]."\n";
				$h_target_drug{lc($drug)}{"flag"} = 1;
				foreach my $var( split(/;;/,$var_list) ){
					$h_all_variants{$var}++;
				}
			}
		}
	}
}
close IN;

open OUT,">$ARGV[1].detail.xls" or die $!;
foreach my $drug( sort {$a cmp $b} keys(%h_target_drug)){
	
	if( $h_target_drug{$drug}{"flag"} == 0){
		print STDERR "Err: no info for drug $drug\n";
		next;
	}
	my $str_drug = $h_target_drug{$drug}{"info"};
	my $lev_1_2 = 0;
	my $lev_3 = 0;
	foreach my $var( sort {$a cmp $b} keys(%{$h_target_drug{$drug}{"var"}})){
		#print "\t$var";
		my $str_var = $var;
		foreach my $allele( sort {$a cmp $b} keys(%{$h_target_drug{$drug}{"var"}{$var}})){

			#只保留高等级的变异
			$lev_1_2 = 1 if($var =~ /^[12]/);
			if($var =~ /^3/){
				$lev_3 = 1;
				next if( $lev_1_2 == 1 );
			}
			if($var =~ /^4/){
				next if( $lev_1_2 == 1 || $lev_3 == 1 );
			}

			print OUT "$str_drug\t$str_var\t$allele\n";
			$str_var = "\t\t\t\t\t";
			$str_drug = "";
		}
		$str_drug = "";
	}
}
close OUT;
open OUT,">$ARGV[1].var.xls" or die $!;
foreach my $var( sort {$a cmp $b} keys(%h_all_variants)){

	if(defined($h_hap{$var})){
		foreach my $mut( sort {$a cmp $b} keys(%{$h_hap{$var}})){
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

