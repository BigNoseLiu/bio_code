#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use MongoDB;
use Encode;
use BSON;
#use BSON::Doc;
use BSON::Types ':all';


my $db_col_name = $ARGV[0];
my $task_id = $ARGV[1];


   #host => '1.12.236.215:27017',
my $client = MongoDB::MongoClient->new(
   host => 'gmbzero.tpddns.cn:31112',
   username => 'lintop',
   password => 'lintop.hx321.mongo'
);
$client->connect;


my $collection = $client->ns( $db_col_name );



my $line_header = <STDIN>;
chomp $line_header;
my %h_header = ();
my @arr_header = split(/\t/,$line_header);
my %h_count = ();
for(my$i=0;$i<@arr_header;$i++){
	my $t_header = $arr_header[$i];
	$t_header =~ s/_SCORE/_score/;
	$t_header = "GERP++_score" if($t_header eq "GERP++_RS");
	$t_header = "CADD_score" if($t_header eq "CADD_raw");
	$t_header = "Eigen_score" if($t_header eq "Eigen-PC-raw_coding");
	if( defined($h_header{$t_header}) ){
		$h_count{$t_header}++;
		$t_header .= "_".$h_count{$t_header};
	}
	$h_header{$t_header} = $i;
	#print STDERR $t_header."\n";
	$h_count{$t_header}++;
}


print `date`;
my $line;
my %h_temp = ();
my %h_mongo_title = ( "out"=>{'AF_popmax'=>'1','AF_male'=>'2'});
my %h_in = ( 'task_id' => $task_id );
my @muts = ();
while( $line = <STDIN> ){
	chomp $line;
	my @arr = split(/\t/,$line);
	my %h_record = ('task_id'=>$task_id);
	#sample info
	$h_record{'sample'}{'id'} = "S100011";

	#init tags for filter
	$h_record{'tags'}{'reported_path'} = 0;
	$h_record{'tags'}{'freq_1'} = 0;
	$h_record{'tags'}{'freq_5'} = 0;
	$h_record{'tags'}{'strong_mut'} = 0;
	#$h_record{'tags'}{'related_path'} = 0;
	#$h_record{'tags'}{'second_finding'} = 0;
	#$h_record{'tags'}{'report_path'} = 0;
	#$h_record{'tags'}{'denovo_mut'} = 0;
	#$h_record{'tags'}{'match_mode'} = 0;
	
	#location
	$h_record{'location'}{'hg19'}{'chr'} = $arr[$h_header{'Chr'}];
	$h_record{'location'}{'hg19'}{'pos1'} = $arr[$h_header{'Start'}];
	$h_record{'location'}{'hg19'}{'pos2'} = $arr[$h_header{'End'}];
	$h_record{'location'}{'hg19'}{'ref'} = $arr[$h_header{'Ref'}];
	$h_record{'location'}{'hg19'}{'alt'} = $arr[$h_header{'Alt'}];



	#result gatk
	my %h_result = ();
	my $gatk_info = $arr[-3];
	my $standard_result = '';
	$h_result{'result_en'} = 'wild';
	$h_result{'result'} = decode_utf8('野生型');
	if(  $gatk_info !~ /AF=/ || $gatk_info =~ /AF=0\.5\d+;/ ){
		$h_result{'result_en'} = 'het';
		$h_result{'result'} = decode_utf8('杂合');
	}
	elsif( $gatk_info =~ /AF=1\.\d+;/ ){
		$h_result{'result_en'} = 'hom';
		$h_result{'result'} = decode_utf8('纯合');
	}
	else{
		print STDERR "Error: illegal AF to define result: $line\n";
	}
	$h_record{'test_result'}{'standard_result'} = $h_result{'result'};
	$h_record{'test_result'}{'standard_result_en'} = $h_result{'result_en'};

	$h_result{'detail'} = $gatk_info;
	$h_result{'detail'} =~ s/^(.+;)ANN_STANDARD=.*$/$1/;
	$h_result{'info'} = "$arr[-4]|$arr[-1]";
	my %h_sample_result = ();
	my @arr_result = ();
	$h_sample_result{'result1'} = \%h_result;
	$h_sample_result{'result2'} = \%h_result;
	$h_sample_result{'sample_id'} = 's2222332';
	push @arr_result,\%h_sample_result;
	$h_sample_result{'sample_id'} = 'S18880';
	push @arr_result,\%h_sample_result;
	
	$h_record{'test_result'}{'mut_source'} = decode_utf8('未定义');
	$h_record{'test_result'}{'note'} = decode_utf8("PS2:患者的新发变异，且无家族史(经双亲验证).\nPM3:在隐性遗传病中，在反式位置上检测到致病变异.\nPM6;未经父母样本验证的新发变异.\nPP1:突变与疾病在家系中共分离(在家系多个患者中检测到此变异)）注:如有更多的证据，可作为更强的证据.\nBS4:在一个家系成员中缺乏共分离.");
	$h_record{'test_result'}{'results'} = \@arr_result;

	#dbsnp
	$h_record{'dbsnp'} = $arr[$h_header{'avsnp150'}];


	my %h_known_var = ();
	#HGMD	POS="1:865595:865595";HGMD_ID="CM1613956";Mutation="SAMD11:NM_152486:c.133A>G:p.K45E";Pathogenicity="DM?";Disease="Retinitis_pigmentosa";Pubmed="27734943"
	#Clinvar	
	$h_known_var{'same_pos'}{'hgmd'}{$arr[$h_header{'HGMD2017_samepos'}]}++;
	$h_known_var{'flank'}{'hgmd'}{$arr[$h_header{'HGMD2017_flank5'}]}++;
	$h_known_var{'flank'}{'clinvar'}{$arr[$h_header{'clinvar_20220416_flank5'}]}++;
	$h_known_var{'same_mut'}{'clinvar'}{$arr[$h_header{'clinvar_anno'}]}++;

	foreach my $known_var_type( sort {$a cmp $b} keys(%h_known_var) ){
		#if(defined($h_known_var{$known_var_type}{'hgmd'})){
		foreach my $db_name(sort {$a cmp $b} keys(%{$h_known_var{$known_var_type}})){
			foreach my $multi_anno(sort {$a cmp $b} keys(%{$h_known_var{$known_var_type}{$db_name}})){
				$multi_anno =~ s/,POS=/&&&,POS=/g;
				foreach my $var_anno(split(/&&&/,$multi_anno)){
					$var_anno .=";";
					if( $var_anno =~ /POS="?([^;"]+)"?;/ ){
						my $temp_pos = $1;
						my $mut_id = $temp_pos;
						$mut_id =~ s/:/_/g;
						$mut_id =~ s/-/X/g;
						$mut_id =~ s/\./X/g;
						$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'pos'} = $temp_pos;
						if( $db_name eq "hgmd" ){
							if( $var_anno =~ /Mutation="([^;]+)";/ ){
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'mut_name'} = $1;
							}
							if( $var_anno =~ /Pathogenicity="([^;]+)";/ ){
								my $t_hgmd_path = $1;
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'path'} = $t_hgmd_path;
							}
							if( $var_anno =~ /Disease="([^;]+)";/ ){
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'disease'} = $1;
							}
							if( $var_anno =~ /Pubmed="([^;]+)";/ ){
								my $pub_id = $1;
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'ref'}{"pubmed_$pub_id"}{'type'} = "pubmed";
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'ref'}{"pubmed_$pub_id"}{'name'} = "pubmed:$pub_id";
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'ref'}{"pubmed_$pub_id"}{'link'} = "https://pubmed.ncbi.nlm.nih.gov/$pub_id";
							}
						}
						if( $db_name eq "clinvar" ){
							if( $var_anno =~ /CLNHGVS=([^;]+);/ ){
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'mut_name'} = $1;
							}
							if( $var_anno =~ /CLNSIG=([^;]+);/ ){
								my $t_clinvar_path = $1;
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'path'} = $t_clinvar_path;
								$h_record{'known_vars'}{'report_path_stat'}{'name'} = $t_clinvar_path if($known_var_type eq "same_mut");
								#$h_record{'tags'}{'reported_path'} = 1 if(defined($h_record{'known_vars'}{'report_path_stat'}{'name'}) && $h_record{'known_vars'}{'report_path_stat'}{'name'} =~ /path/i);
							}
							if( $var_anno =~ /CLNREVSTAT=([^;]+);/ ){
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'path_level'} = $1;
							}
							if( $var_anno =~ /CLNDISDB=([^;]+);/ ){
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'disease'} = $1;
							}
							if( $var_anno =~ /CLNVARID=([^;]+);/ ){
								my $clinvar_id = $1;
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'ref'}{"clinvarid_$clinvar_id"}{'type'} = "clinvar";
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'ref'}{"clinvarid_$clinvar_id"}{'name'} = "ClinVar";
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'ref'}{"clinvarid_$clinvar_id"}{'link'} = "https://www.ncbi.nlm.nih.gov/clinvar/variation/$clinvar_id";
							}
						}
					}
				}
			}
		}
	}
	#only the same
	$h_record{'tags'}{'reported_path'} = 1 if(defined($h_record{'known_vars'}{'report_path_stat'}{'name'}) && $h_record{'known_vars'}{'report_path_stat'}{'name'} =~ /path/i);


	#hgvs
	my $t_line = $line." ";
        if( $t_line =~ /ANN=([^;\t]+)[;\s]/ ){
                my $t_snpeff_annos = $1;
		#get standrad transc
		my $standard_trans = "-";
        	if( $t_line =~ /ANN_STANDARD=([^;\t]+)[;\s]/ ){
                        my @t_arr = split(/\|/,$1);
			$standard_trans = $t_arr[6];
		}
		
		my @arr_hgvs = ();
                foreach my $t_snpeff_anno( split(/,/, $t_snpeff_annos) ){
                        my @t_arr = split(/\|/,$t_snpeff_anno);
                        my $t_trans_type = "-";
                        $t_arr[2] = "-" if( !defined($t_arr[2]) );
                        $t_arr[6] = "-" if( !defined($t_arr[6]) );
			$h_record{'tags'}{'strong_mut'} = 1 if($t_arr[2] =~ /HIGH/i);

			my %h_hgvs = ();
			$h_hgvs{'gene'} = $t_arr[3];
			#get the genes
			#$h_record{'genes'}{$t_arr[3]}{'flag'} = 1 if($t_arr[3] =~ /\S/);
			$h_hgvs{'trans'} = $t_arr[6];
			$h_hgvs{'standard'} = 0;
			if( $t_arr[6] eq $standard_trans ){
				$h_hgvs{'standard'} = 1;
			}
			if( defined($t_arr[11]) ){
				if( $t_arr[11] =~ /^\d+\/(\d+)$/ ){
					$h_hgvs{'len'} = $1;
				}
			}
			if( defined($t_arr[8]) ){
				if( $t_arr[8] =~ /^(\d+)\/(\d+)$/ ){
					$h_hgvs{'exon_intron'}{'locate'} = $1;
					$h_hgvs{'exon_intron'}{'total'} = $2;
					if( $t_arr[1] =~ /intron/i || $t_arr[1] =~ /splice_donor_variant/i ){
						$h_hgvs{'exon_intron'}{'type'} = 'I';
					}
					else{
						$h_hgvs{'exon_intron'}{'type'} = 'E';
					}
				}
			}
			$h_hgvs{'c_point'} = $t_arr[9] if(defined($t_arr[9]) && $t_arr[9] =~ /\S/);
			$h_hgvs{'p_point'} = $t_arr[10] if(defined($t_arr[10]) && $t_arr[10] =~ /\S/);
			if( $t_arr[1] =~ /\S/ ){
				my $t_count = 0;
				foreach my $t_type(split(/&/,$t_arr[1])){
					$t_count ++;
					$h_hgvs{'type'}{$t_type}{'ch_name'} = $t_type;
					$h_hgvs{'type'}{$t_type}{'priority'} = $t_count;
				}
			}
			push @arr_hgvs,\%h_hgvs;

                }
		$h_record{'hgvs'} = \@arr_hgvs;
	}



=pod
	my @arr_disease = ();
	my %h_disease = ();
	$h_disease{'id'} = 'DI10123';
	$h_disease{'omim_id'} = '101123';
	$h_disease{'orphanet_id'} = 'ORP11123';
	$h_disease{'en_name'} = 'hearing loss type 1A';
	$h_disease{'ch_name'} = decode_utf8('遗传性耳聋1A');
	$h_disease{'morbidity'} = decode_utf8('1/10万');
	$h_disease{'age_onset'} = decode_utf8('青少年');
	$h_disease{'inher_mode'}{'AR'}{'total'} = decode_utf8('常染色体隐性');
	$h_disease{'inher_mode'}{'AR'}{'short'} = decode_utf8('常隐');
	$h_disease{'inher_mode'}{'AD'}{'total'} = decode_utf8('常染色体显性');
	$h_disease{'inher_mode'}{'AD'}{'short'} = decode_utf8('常显');
	$h_disease{'descr'} = decode_utf8('青少年');
	push @arr_disease,\%h_disease;
	push @arr_disease,\%h_disease;
	$h_record{'disease'}{'BRCA1'} = \@arr_disease;
	$h_record{'disease'}{'GJB2'} = \@arr_disease;

	#clinvar stat
	#my @arr_clinvarstat= ();
	#my %h_clinvarstat= ();
	$h_record{'clinvar_stat'}{'BRCA1'}{'nonsense'}{'path'} = 10;
	$h_record{'clinvar_stat'}{'BRCA1'}{'nonsense'}{'likelypath'} = 21;
	$h_record{'clinvar_stat'}{'BRCA1'}{'nonsense'}{'vus'} = 1;
	$h_record{'clinvar_stat'}{'BRCA1'}{'nonsense'}{'likelybenign'} = 3;
	$h_record{'clinvar_stat'}{'BRCA1'}{'nonsense'}{'benign'} = 7;
	$h_record{'clinvar_stat'}{'BRCA1'}{'missense'}{'benign'} = 7;
	$h_record{'clinvar_stat'}{'GJB2'}{'nonsense'}{'benign'} = 7;
=cut

	#phenotype
	$h_record{'phenotype'}{'note'} = decode_utf8('PP4:变异携带者的表型或家族史高度符合某种单基因遗传病');

	my %h_pheno = ();
	$h_pheno{'source'} = 'HPO';
	$h_pheno{'hpo_id'} = 'HP:101011';
	$h_pheno{'match_flag'} = 1;
	$h_pheno{'match_type'} = '精准';
	$h_pheno{'gene_uniq'} = 'uniq';
	$h_pheno{'ch_name'} = decode_utf8('耳聋1');
	$h_pheno{'en_name'} = 'deafness,hearing impair 1';
	$h_record{'phenotype'}{'ID_001101'} = \%h_pheno;

	#popfreq
	$h_record{'popfreq'}{'note'} = decode_utf8("待修改");
	my $freq_database = 'gnomad_genome';
	$h_record{'popfreq'}{'detail'}{$freq_database}{'all'}{'freq'} = $arr[$h_header{'AF'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'afr'}{'freq'} = $arr[$h_header{'AF_afr'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'sas'}{'freq'} = $arr[$h_header{'AF_sas'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'amr'}{'freq'} = $arr[$h_header{'AF_amr'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'eas'}{'freq'} = $arr[$h_header{'AF_eas'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'nfe'}{'freq'} = $arr[$h_header{'AF_nfe'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'fin'}{'freq'} = $arr[$h_header{'AF_fin'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'asj'}{'freq'} = $arr[$h_header{'AF_asj'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'oth'}{'freq'} = $arr[$h_header{'AF_oth'}];
	$freq_database = 'gnomad_exome';
	$h_record{'popfreq'}{'detail'}{$freq_database}{'all'}{'freq'} = $arr[$h_header{'AF_2'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'afr'}{'freq'} = $arr[$h_header{'AF_afr_2'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'sas'}{'freq'} = $arr[$h_header{'AF_sas_2'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'amr'}{'freq'} = $arr[$h_header{'AF_amr_2'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'eas'}{'freq'} = $arr[$h_header{'AF_eas_2'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'nfe'}{'freq'} = $arr[$h_header{'AF_nfe_2'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'fin'}{'freq'} = $arr[$h_header{'AF_fin_2'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'asj'}{'freq'} = $arr[$h_header{'AF_asj_2'}];
	$h_record{'popfreq'}{'detail'}{$freq_database}{'oth'}{'freq'} = $arr[$h_header{'AF_oth_2'}];
	
	$h_record{'popfreq'}{'used_freq'}{"freq"} = -1;
	$h_record{'popfreq'}{'used_freq'}{"source"} = "not_detected";
	$h_record{'popfreq'}{'used_freq'}{"rule"} = "max";

	foreach my $t_db( keys(%{ $h_record{'popfreq'}{'detail'} }) ){
		foreach my $t_pop( keys(%{ $h_record{'popfreq'}{'detail'}{$t_db} }) ){
			#change . to -1 for population freq
			if( $h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'} eq '.' ){
				$h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'} = -1;
			}
			#get the max freq
			if($h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'} > $h_record{'popfreq'}{'used_freq'}{"freq"}){
				$h_record{'popfreq'}{'used_freq'}{"freq"} = $h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'};
				$h_record{'popfreq'}{'used_freq'}{"source"} = "$t_db:$t_pop";
			}
			#reserve for hom/het count
			$h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'hom'} = -1;
			$h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'het'} = -1;
		}
	}
	$h_record{'tags'}{'freq_1'} = 1 if( $h_record{'popfreq'}{'used_freq'}{"freq"} <=0.01 );
	$h_record{'tags'}{'freq_5'} = 1 if( $h_record{'popfreq'}{'used_freq'}{"freq"} <=0.05 );

	#预测
	#$h_record{'predicted'}{'note'} = encode("utf-8",decode("GB2312","PP3:多种统计方法预测出该变异会对基因或基因产物造成有害的影响，包括保守性预测、进化预测、剪接位点影响等\nBP4:多种统计方法预测出该变异会对基因或基因产物无影响，包括保守性预测、进化预测、剪接位点影响等。\nBP7:同义变异且预测不影响剪接。"));
	$h_record{'predicted'}{'note'} = decode_utf8("PP3:多种统计方法预测出该变异会对基因或基因产物造成有害的影响，包括保守性预测、进化预测、剪接位点影响等\nBP4:多种统计方法预测出该变异会对基因或基因产物无影响，包括保守性预测、进化预测、剪接位点影响等。\nBP7:同义变异且预测不影响剪接。");
	my @arr_predict = ('REVEL','VEST4','SIFT','SIFT4G','Polyphen2_HDIV','Polyphen2_HVAR','LRT','MutationTaster','FATHMM','PROVEAN','MetaSVM','MetaLR','MetaRNN','M-CAP','PrimateAI','DEOGEN2');
	my %h_predict = ();
	$h_predict{'综合评估'}{'BayesDel_addAF'}{'note'} = '依据：多个预测工具和人群频率';
	$h_predict{'综合评估'}{'BayesDel_noAF'}{'note'} = '依据：多个预测工具和人群频率';
	$h_predict{'综合评估'}{'BayesDel_noAF'}{'note'} = '依据：多个预测工具和人群频率';
	$h_predict{'综合评估'}{'ClinPred'}{'note'} = '依据：多个预测工具和人群频率';
	$h_predict{'综合评估'}{'CADD'}{'note'} = '依据：模拟突变自然选择规律';
	$h_predict{'综合评估'}{'DANN'}{'note'} = '依据：模拟突变自然选择规律';
	$h_predict{'综合评估'}{'DEOGEN2'}{'note'} = '依据：预测工具及其他多种特征';
	$h_predict{'综合评估'}{'M-CAP'}{'note'} = '依据：多个预测工具和物种序列比对';
	$h_predict{'综合评估'}{'MetaLR'}{'note'} = '依据：多个预测工具和人群频率';
	$h_predict{'综合评估'}{'MetaRNN'}{'note'} = '依据：多个预测工具';
	$h_predict{'综合评估'}{'MetaSVM'}{'note'} = '依据：多个预测工具和人群频率';
	$h_predict{'综合评估'}{'REVEL'}{'note'} = '依据：多个预测工具';
	$h_predict{'综合评估'}{'Eigen'}{'note'} = '依据：多个预测工具和人群频率';
	$h_predict{'保守性'}{'GERP++'}{'note'} = '依据：基因组进化速率分析';
	$h_predict{'基因组功能'}{'GenoCanyon'}{'note'} = '依据：多个预测工具和实验数据';
	#$h_predict{'线粒体'}{'MitoTip'}{'note'} = '依据：人群频率、核苷酸变化性质和保护';
	#$h_predict{'剪接'}{'SpliceAI'}{'note'} = '依据：深度神经网络';
	$h_predict{'剪接'}{'dbscSNV_ADA'}{'note'} = '';
	$h_predict{'剪接'}{'dbscSNV_RF'}{'note'} = '';
	$h_predict{'编码功能'}{'FATHMM'}{'note'} = '依据：进化保守性';
	$h_predict{'编码功能'}{'LIST-S2'}{'note'} = '依据：进化保守性';
	$h_predict{'编码功能'}{'LRT'}{'note'} = '依据：进化保守性';
	$h_predict{'编码功能'}{'MPC'}{'note'} = '依据：所在位置';
	$h_predict{'编码功能'}{'MVP'}{'note'} = '依据：机器学习';
	$h_predict{'编码功能'}{'MutPred'}{'note'} = '依据：预测工具及其他多种特征';
	$h_predict{'编码功能'}{'MutationAssessor'}{'note'} = '依据：进化保守性';
	$h_predict{'编码功能'}{'MutationTaster'}{'note'} = '依据：蛋白结构/功能和进化保守性';
	$h_predict{'编码功能'}{'PROVEAN'}{'note'} = '依据：蛋白序列比对';
	$h_predict{'编码功能'}{'Polyphen2_HDIV'}{'note'} = '依据：蛋白结构/功能和进化保守性(复杂疾病相关)';
	$h_predict{'编码功能'}{'Polyphen2_HVAR'}{'note'} = '依据：蛋白结构/功能和进化保守性(孟德尔遗传病相关)';
	$h_predict{'编码功能'}{'PrimateAI'}{'note'} = '依据：进化保守性';
	$h_predict{'编码功能'}{'SIFT'}{'note'} = '依据：进化保守性';
	$h_predict{'编码功能'}{'SIFT4G'}{'note'} = '依据：进化保守性';
	$h_predict{'编码功能'}{'VEST4'}{'note'} = '依据：机器学习';
	$h_predict{'编码功能'}{'fathmm-MKL_coding'}{'note'} = '依据：进化保守性（预测编码及非编码SNV）';
	$h_predict{'编码功能'}{'fathmm-XF_coding'}{'note'} = '依据：进化保守性（预测编码及非编码SNV）';
	$h_predict{'编码功能'}{'integrated_fitCons'}{'note'} = '依据：基于进化的潜在基因组功能';

	#'REVEL'=>"",'VEST4','SIFT','SIFT4G','Polyphen2_HDIV','Polyphen2_HVAR','LRT','MutationTaster','FATHMM','PROVEAN','MetaSVM','MetaLR','MetaRNN','M-CAP','PrimateAI','DEOGEN2');
	$h_record{'predicted'}{'damage_count'} = 0;
	$h_record{'predicted'}{'predicted_count'} = 0;
	foreach my $class(keys(%h_predict)){
		#foreach my $t_db(@arr_predict){
		foreach my $t_db(keys(%{$h_predict{$class}})){
			if( defined($h_header{$t_db.'_score'}) ){
				$h_record{'predicted'}{'total_count'}++;
				$arr[$h_header{$t_db.'_score'}] = -1 if($arr[$h_header{$t_db.'_score'}] eq '.');
				$h_record{'predicted'}{'detail'}{$t_db}{'type'} = decode_utf8($class);
				$h_record{'predicted'}{'detail'}{$t_db}{'order'} = 1;
				$h_record{'predicted'}{'detail'}{$t_db}{'note'} = decode_utf8($h_predict{$class}{$t_db}{'note'});
				$h_record{'predicted'}{'detail'}{$t_db}{'score'} = $arr[$h_header{$t_db.'_score'}];
				if( defined($h_header{$t_db.'_pred'}) ){
					my $t_pred = $arr[$h_header{$t_db.'_pred'}];
					#print STDERR $t_pred."\n";
					$h_record{'predicted'}{'detail'}{$t_db}{'pred'} = $t_pred;
			
					if( $t_pred ne '.' ){
						$h_record{'predicted'}{'predicted_count'}++;
					}
					if( $t_pred =~ /D/ ){
						$h_record{'predicted'}{'damage_count'}++;
					}
					#$h_temp{$t_db."\t".$arr[$h_header{$t_db.'_pred'}]}++;
				}
			}
			else{
				print STDERR "waring: no predict database $t_db\n";
			}
		}
	}
	$collection->insert_one( \%h_record );
}

$client->disconnect;

#print join("\n",(sort {$a cmp $b} keys(%h_temp)))."\n";
#my $codec = BSON->new;
#my $doc = {'abc'=>"PP3:多种统计方法预测出该变异会对基因或基因产物造成有害的影响，包括保守性预测、进化预测、剪接位点影响等\nBP4:多种统计方法预测出该变异会对基因或基因产物无影响，包括保守性预测、进化预测、剪接位点影>响等。\nBP7:同义变异且预测不影响剪接。", };

#$collection->insert_one( $codec->encode_one($doc) );
#print `date`;
#my @data = $collection->find({"task_id"=>'A10010101'})->all();
#print `date`;
#$h_in{'all_muts'} = \@muts;
#$collection->insert_one( \%h_in);

#0	Chr	:	chr1
#1	Start	:	13116
#2	End	:	13118
#3	Ref	:	TGA
#4	Alt	:	GGG
#5	Func.refGeneWithVer	:	ncRNA_intronic
#6	Gene.refGeneWithVer	:	DDX11L1;LOC102725121
#7	GeneDetail.refGeneWithVer	:	.
#8	ExonicFunc.refGeneWithVer	:	.
#9	AAChange.refGeneWithVer	:	.
#10	CLNALLELEID	:	.
#11	CLNDN	:	.
#12	CLNDISDB	:	.
#13	CLNREVSTAT	:	.
#14	CLNSIG	:	.
#15	avsnp150	:	.
#16	AF	:	.
#17	AF_popmax	:	.
#18	AF_male	:	.
#19	AF_female	:	.
#20	AF_raw	:	.
#21	AF_afr	:	.
#22	AF_sas	:	.
#23	AF_amr	:	.
#24	AF_eas	:	.
#25	AF_nfe	:	.
#26	AF_fin	:	.
#27	AF_asj	:	.
#28	AF_oth	:	.
#29	non_topmed_AF_popmax	:	.
#30	non_neuro_AF_popmax	:	.
#31	non_cancer_AF_popmax	:	.
#32	controls_AF_popmax	:	.
#33	AF	:	.
#34	AF_popmax	:	.
#35	AF_male	:	.
#36	AF_female	:	.
#37	AF_raw	:	.
#38	AF_afr	:	.
#39	AF_sas	:	.
#40	AF_amr	:	.
#41	AF_eas	:	.
#42	AF_nfe	:	.
#43	AF_fin	:	.
#44	AF_asj	:	.
#45	AF_oth	:	.
#46	non_topmed_AF_popmax	:	.
#47	non_neuro_AF_popmax	:	.
#48	non_cancer_AF_popmax	:	.
#49	controls_AF_popmax	:	.
#50	SIFT_score	:	.
#51	SIFT_converted_rankscore	:	.
#52	SIFT_pred	:	.
#53	SIFT4G_score	:	.
#54	SIFT4G_converted_rankscore	:	.
#55	SIFT4G_pred	:	.
#56	LRT_score	:	.
#57	LRT_converted_rankscore	:	.
#58	LRT_pred	:	.
#59	MutationTaster_score	:	.
#60	MutationTaster_converted_rankscore	:	.
#61	MutationTaster_pred	:	.
#62	MutationAssessor_score	:	.
#63	MutationAssessor_rankscore	:	.
#64	MutationAssessor_pred	:	.
#65	FATHMM_score	:	.
#66	FATHMM_converted_rankscore	:	.
#67	FATHMM_pred	:	.
#68	PROVEAN_score	:	.
#69	PROVEAN_converted_rankscore	:	.
#70	PROVEAN_pred	:	.
#71	MetaSVM_score	:	.
#72	MetaSVM_rankscore	:	.
#73	MetaSVM_pred	:	.
#74	MetaLR_score	:	.
#75	MetaLR_rankscore	:	.
#76	MetaLR_pred	:	.
#77	MetaRNN_score	:	.
#78	MetaRNN_rankscore	:	.
#79	MetaRNN_pred	:	.
#80	M-CAP_score	:	.
#81	M-CAP_rankscore	:	.
#82	M-CAP_pred	:	.
#83	MutPred_score	:	.
#84	MutPred_rankscore	:	.
#85	MVP_score	:	.
#86	MVP_rankscore	:	.
#87	MPC_score	:	.
#88	MPC_rankscore	:	.
#89	PrimateAI_score	:	.
#90	PrimateAI_rankscore	:	.
#91	PrimateAI_pred	:	.
#92	DEOGEN2_score	:	.
#93	DEOGEN2_rankscore	:	.
#94	DEOGEN2_pred	:	.
#95	BayesDel_addAF_score	:	.
#96	BayesDel_addAF_rankscore	:	.
#97	BayesDel_addAF_pred	:	.
#98	BayesDel_noAF_score	:	.
#99	BayesDel_noAF_rankscore	:	.
#100	BayesDel_noAF_pred	:	.
#101	ClinPred_score	:	.
#102	ClinPred_rankscore	:	.
#103	ClinPred_pred	:	.
#104	LIST-S2_score	:	.
#105	LIST-S2_rankscore	:	.
#106	LIST-S2_pred	:	.
#107	Aloft_pred	:	.
#108	Aloft_Confidence	:	.
#109	DANN_score	:	.
#110	DANN_rankscore	:	.
#111	fathmm-MKL_coding_score	:	.
#112	fathmm-MKL_coding_rankscore	:	.
#113	fathmm-MKL_coding_pred	:	.
#114	fathmm-XF_coding_score	:	.
#115	fathmm-XF_coding_rankscore	:	.
#116	fathmm-XF_coding_pred	:	.
#117	Eigen-raw_coding	:	.
#118	Eigen-raw_coding_rankscore	:	.
#119	Eigen-PC-raw_coding	:	.
#120	Eigen-PC-raw_coding_rankscore	:	.
#121	integrated_fitCons_score	:	.
#122	integrated_fitCons_rankscore	:	.
#123	integrated_confidence_value	:	.
#124	GERP++_NR	:	.
#125	GERP++_RS	:	.
#126	GERP++_RS_rankscore	:	.
#127	phyloP100way_vertebrate	:	.
#128	phyloP100way_vertebrate_rankscore	:	.
#129	phyloP30way_mammalian	:	.
#130	phyloP30way_mammalian_rankscore	:	.
#131	phastCons100way_vertebrate	:	.
#132	phastCons100way_vertebrate_rankscore	:	.
#133	phastCons30way_mammalian	:	.
#134	phastCons30way_mammalian_rankscore	:	.
#135	SiPhy_29way_logOdds	:	.
#136	SiPhy_29way_logOdds_rankscore	:	.
#137	Interpro_domain	:	.
#138	GTEx_V8_gene	:	.
#139	GTEx_V8_tissue	:	.
#140	dbscSNV_ADA_SCORE	:	.
#141	dbscSNV_RF_SCORE	:	.
#142	Otherinfo1	:	0.5
#143	Otherinfo2	:	104.64
#144	Otherinfo3	:	24
#145	Otherinfo4	:	chr1
#146	Otherinfo5	:	13116
#147	Otherinfo6	:	.
#148	Otherinfo7	:	TGA
#149	Otherinfo8	:	GGG
#150	Otherinfo9	:	104.64
#151	Otherinfo10	:	.
#152	Otherinfo11	:	AC=1;AF=0.500;AN=2;BaseQRankSum=0.159;ClippingRankSum=0.00;DP=26;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=38.73;MQRankSum=-2.061e+00;QD=4.36;RAW_MQ=39000.00;ReadPosRankSum=-3.124e+00;SOR=0.368
#153	Otherinfo12	:	GT:AD:DP:GQ:PL
#154	Otherinfo13	:	0/1:19,5:24:99:112,0,783
