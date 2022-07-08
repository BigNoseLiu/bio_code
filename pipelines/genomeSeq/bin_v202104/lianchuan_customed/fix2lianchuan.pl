my @types = ("other","het","hom");
my $call_type = "somatic";
my $line;
my $line_count = 0;
#my @headers = ('Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'esp6500siv2_all', 'gnomAD_genome_ALL', 'gnomAD_genome_AFR', 'gnomAD_genome_EAS', 'gnomAD_genome_AMR', 'avsnp150', 'nci60', 'cosmic81', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'GERP++_RS', 'phyloP20way_mammalian', 'phyloP100way_vertebrate', 'SiPhy_29way_logOdds', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'Filter', 'All_depth', 'Alt_depth', 'Mut_fre');
my @headers = ('Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'esp6500siv2_all', '1000g2014oct_all', '1000g2014oct_afr', '1000g2014oct_eas', '1000g2014oct_eur', 'avsnp150', 'nci60', 'cosmic81', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'RadialSVM_score', 'RadialSVM_pred', 'LR_score', 'LR_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'GERP++_RS', 'phyloP46way_placental', 'phyloP100way_vertebrate', 'SiPhy_29way_logOdds', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG', 'Filter', 'All_depth', 'Alt_depth', 'Mut_fre');
print join("\t",@headers)."\n";
my %h_header = ();
while( $line = <STDIN> ){
	chomp $line;
	$line_count++;
	my @arr = split(/\t/,$line);
	if( $line_count == 1 ){
		$line =~ s/WithVer//g;
		$line =~ s/cosmic70/cosmic81/;
		my @t_headers = split( /\t/, $line );
		for( $i=0;$i<scalar@t_headers;$i++ ){
			if( $t_headers[$i] eq "Otherinfo10"){
				$t_headers[$i] = "Filter";
			}
			$h_header{$t_headers[$i]} = $i;
		}
		next;
	}
	#next if( $arr[$h_header{"Filter"}] !~/PASS/i );
	$arr[$h_header{"Filter"}] = "PASS" if($line_count > 1); 
	$arr[0] = "chr".$arr[0];
	for( $i=0;$i<(scalar@headers-3);$i++ ){
		if( !defined($h_header{$headers[$i]}) ){
			print "-\t";
		}
		else{
			print $arr[$h_header{$headers[$i]}]."\t";
		}
	}
	
	my $type = $types[0];
	if( $arr[$h_header{"Otherinfo13"}] =~ /^0\/1:/i || $arr[$h_header{"Otherinfo13"}] =~ /^1\/0:/i ){
		$type = $types[1];
	}
	elsif( $arr[$h_header{"Otherinfo13"}] =~ /^1\/1:/i ){
		$type = $types[2];
	}
	if( $call_type =~ /somatic/i ){
		my %h_format_ids = ();
		my $format_id_count = 0;
		foreach my $format_id( split(/:/,$arr[$h_header{"Otherinfo12"}]) ){
			$h_format_ids{$format_id} = $format_id_count;
			$format_id_count++;
		}
		my @datas =  split(/:/,$arr[$h_header{"Otherinfo13"}]);
		my $dp = $datas[$h_format_ids{"DP"}];
		my @ads = split(/,/,$datas[$h_format_ids{"AD"}]);
		my $ad = $ads[-1];
		my $af = ($ad/$dp)*100;
		$af =~ s/(\.\d\d)\d*$/$1/;
		$type = $af;
		print "$dp\t$ad\t$af%\n";
	}
}
