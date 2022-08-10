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
my $main_sample = $ARGV[2];


   #host => '1.12.236.215:27017',
   #host => 'gmbzero.tpddns.cn:31112',
   #host => '10.10.9.35:31112',
my $client = MongoDB::MongoClient->new(
   host => '10.10.9.35:31112',
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
	$h_header{$t_header} = $i;
}


my %h_mut_type2name = ('missense_variant'=> decode_utf8('错义'),'synonymous_variant'=> decode_utf8('同义'),'frameshift_variant'=> decode_utf8('移码'),'intron_variant'=> decode_utf8('内含子'),'splice_region_variant'=> decode_utf8('剪切'),'3_prime_UTR_variant'=> decode_utf8('UTR3'),'5_prime_UTR_variant'=> decode_utf8('UTR5'),'downstream_gene_variant'=> decode_utf8('下游'),'upstream_gene_variant'=> decode_utf8('上游'),'splice_region_variant'=> decode_utf8('剪切区域'),'intergenic_region'=> decode_utf8('基因间区'), 'non_coding_transcript_exon_variant'=> decode_utf8('外显子'), 'intragenic_variant'=>decode_utf8('基因区'));
print "start at:\t".`date`;
my $line;
my %h_temp = ();
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
	my ( $t_chr,$t_pos1,$t_pos2,$t_ref,$t_alt ) = ($arr[$h_header{'Chr'}],$arr[$h_header{'Start'}],$arr[$h_header{'End'}],$arr[$h_header{'Ref'}],$arr[$h_header{'Alt'}]);
	$h_record{'location'}{'hg19'}{'chr'} = $t_chr;
	$h_record{'location'}{'hg19'}{'pos1'} = $t_pos1;
	$h_record{'location'}{'hg19'}{'pos2'} = $t_pos2;
	$h_record{'location'}{'hg19'}{'ref'} = $t_ref;
	$h_record{'location'}{'hg19'}{'alt'} = $t_alt;
	if( $t_chr =~ /^chrm/i || $t_chr =~ /^m/i ){
		$h_record{'location'}{'hg19'}{'name'} = "m.$t_pos1$t_ref>$t_alt";
		if( $t_ref=~/^[ATCG]$/i && $t_alt=~/^[ATCG]$/i ){
		}
		elsif( $t_alt eq "-" ){
			if(length($t_ref)==1){
				$h_record{'location'}{'hg19'}{'name'} = "m.".$t_pos1."del";
			}
			else{
				$h_record{'location'}{'hg19'}{'name'} = "m.$t_pos1".'_'.($t_pos1+length($t_ref)-1)."del";
			}
		}
		elsif( $t_ref eq "-" ){
			$h_record{'location'}{'hg19'}{'name'} = "m.".$t_pos1.'_'.($t_pos1+length($t_ref))."ins".$t_alt;
		}
		else{
			$h_record{'location'}{'hg19'}{'name'} = "m.".$t_pos1.'_'.($t_pos1+length($t_ref))."delins".$t_alt;
		}
	}

	#result gatk
	$h_record{'test_result'}{'standard_result_en'} = "-";
	$h_record{'test_result'}{'standard_result'} = "-";

	#获取各样本数据
	my @arr_result = ();
	foreach my $af_header( keys(%h_header) ){
		my %h_sample_result = ();
		if( $af_header =~ /^mut_result_(\S+)$/ ){
			my $sample_id = $1;
			my %h_result = ();
			$h_result{'result_en'} = '.';
			$h_result{'result'} = ".";
			if(  $arr[$h_header{$af_header}] =~ /het/ ){
				$h_result{'result_en'} = 'het';
				$h_result{'result'} = decode_utf8('杂合');
			}
			elsif(  $arr[$h_header{$af_header}] =~ /hom/ ){
				$h_result{'result_en'} = 'hom';
				$h_result{'result'} = decode_utf8('纯合');
			}
			elsif(  $arr[$h_header{$af_header}] =~ /\d/ ){
				$h_result{'result_en'} = $arr[$h_header{$af_header}];
				$h_result{'result'} = $arr[$h_header{$af_header}];
			}
			if(!defined($main_sample)){
				$h_record{'test_result'}{'standard_result_en'} .= "|".$h_result{'result_en'};
				$h_record{'test_result'}{'standard_result'}    .= "|".$h_result{'result'};
			}
			elsif( $main_sample eq $sample_id ){
				$h_record{'test_result'}{'standard_result_en'} = $h_result{'result_en'};
				$h_record{'test_result'}{'standard_result'}    = $h_result{'result'};
			}
			$h_result{'info'} = "-";
			$h_result{'detail'} = "-";
			$h_result{'detail'} = $arr[$h_header{"Otherinfo11"}] if(defined($h_header{"Otherinfo11"}));
			$h_result{'detail'} =~ s/^(.+;)ANN_STANDARD=.*$/$1/;
			$h_result{'info'} = $arr[$h_header{"mut_detail_$sample_id"}] if(defined($h_header{"mut_detail_$sample_id"}));
			$h_result{'info'} =~ s/^\.\/\.:.*$/\./;

			$h_sample_result{'result1'} = \%h_result;
			$h_sample_result{'result2'} = \%h_result;
			$h_sample_result{'sample_id'} = $sample_id;
			push @arr_result,\%h_sample_result;
		}
	}

	$h_record{'test_result'}{'standard_result_en'} =~ s/^-\|//;
	$h_record{'test_result'}{'standard_result'} =~  s/^-\|//;

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
	$h_known_var{'flank'}{'clinvar'}{$arr[$h_header{'clinvar_flank5_20220416'}]}++;
	$h_known_var{'flank'}{'clinvar'}{$arr[$h_header{'clinvar_samepos_20220416'}]}++;
	$h_known_var{'same_mut'}{'clinvar'}{$arr[$h_header{'clinvar_anno'}]}++;

	my %h_known_path = ();
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
						$h_record{'known_vars'}{$known_var_type}{$mut_id}{'pos'} = $temp_pos;
						$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'pos'} = $temp_pos;
						if( $db_name eq "hgmd" ){
							if( $var_anno =~ /Mutation="([^;]+)";/ ){
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'mut_name'} = $1;
							}
							if( $var_anno =~ /Pathogenicity="([^;]+)";/ ){
								my $t_hgmd_path = $1;
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'path'} = $t_hgmd_path;
								$h_known_path{$known_var_type}{$t_hgmd_path}++;
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
								#$h_record{'known_vars'}{$known_var_type}{$mut_id}{'pos'} = $1;
							}
							if( $var_anno =~ /CLNSIG=([^;]+);/ ){
								my $t_clinvar_path = $1;
								my $t_clinvar_path_raw = $t_clinvar_path;
								if($t_clinvar_path =~ /Likely_benign/i){
									$t_clinvar_path = decode_utf8('可能良性');
								}
								elsif($t_clinvar_path =~ /benign/i){
									$t_clinvar_path = decode_utf8('良性');
								}
								elsif($t_clinvar_path =~ /Uncertain_significance/i){
									$t_clinvar_path = decode_utf8('VUS');
								}
								elsif($t_clinvar_path =~ /Conflicting_interpretations_of_pathogenicity/i){
									$t_clinvar_path = decode_utf8('冲突');
								}
								elsif($t_clinvar_path =~ /Likely_pathogenic/i){
									$t_clinvar_path = decode_utf8('可能致病');
								}
								elsif($t_clinvar_path =~ /Pathogenic/i){
									$t_clinvar_path = decode_utf8('致病');
								}
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'path'} = $t_clinvar_path;
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'path_raw'} = $t_clinvar_path_raw;
								$h_known_path{$known_var_type}{$t_clinvar_path}++;
								$h_record{'known_vars'}{'report_path_stat'}{'name'} = $t_clinvar_path if($known_var_type eq "same_mut");
								#$h_record{'tags'}{'reported_path'} = 1 if(defined($h_record{'known_vars'}{'report_path_stat'}{'name'}) && $h_record{'known_vars'}{'report_path_stat'}{'name'} =~ /path/i);
							}
							if( $var_anno =~ /CLNREVSTAT=([^;]+);/ ){
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'path_level'} = $1;
							}
							if( $var_anno =~ /CLNDISDB=([^;]+);/ ){
								my $clndisdb = $1;
								$clndisdb =~ s/\/\//,/;
								$h_record{'known_vars'}{$known_var_type}{$mut_id}{$db_name}{'disease'} = $clndisdb;
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
			$h_hgvs{'len'} = '0';
			if( defined($t_arr[11]) ){
				if( $t_arr[11] =~ /^[-\d]+\/(\d+)$/ ){
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
			$h_hgvs{'impact'} = $t_arr[2];
			if( $t_arr[1] =~ /\S/ ){
				my $t_count = 0;
				foreach my $t_type(split(/&/,$t_arr[1])){
					$t_count ++;
					my $t_name = $t_type;
					$t_name = $h_mut_type2name{$t_type} if( defined( $h_mut_type2name{$t_type} ) );
					$h_hgvs{'type'}{$t_type}{'ch_name'} = $t_name;
					$h_hgvs{'type'}{$t_type}{'priority'} = $t_count;
				}
			}
			push @arr_hgvs,\%h_hgvs;

                }
		$h_record{'hgvs'} = \@arr_hgvs;
	}

	#Exomiser 打分
	foreach my $af_header( keys(%h_header) ){
		if( $af_header =~ /^EXOMISER_(\S+)$/ ){
			my $e_type = $1;
			$arr[$h_header{$af_header}] =~ s/\.0*$//;
			$h_record{'rank_score'}{$e_type} = $arr[$h_header{$af_header}];
		}
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
	foreach my $af_header( keys(%h_header) ){
		#线粒体人群频率库
		if( $t_chr =~ /^chrm/i || $t_chr =~ /^m/i ){
			#if( $af_header =~ /^GBFreq_(\S+)$/ ){
			if( $af_header =~ /^mitofreq\.(\S+)\.(\S+)$/ ){
				my $t_db_type = $1;
				my $t_af_type = $2;
				$h_record{'popfreq'}{'detail'}{'mitofreq'}{$t_db_type}{$t_af_type} = $arr[$h_header{$af_header}];
			}
		}
		else{#非线粒体人群频率库
			if( $af_header =~ /^ExAC_(\S+)$/ ){
				$h_record{'popfreq'}{'detail'}{'ExAC'}{$1}{'freq'} = $arr[$h_header{$af_header}];
			}
			elsif( $af_header =~ /^1000g2015aug_(\S+)$/ ){
				$h_record{'popfreq'}{'detail'}{'1000genome'}{$1}{'freq'} = $arr[$h_header{$af_header}];
			}
			elsif( $af_header =~ /^gnomadExome_(\S+)$/ ){
				$h_record{'popfreq'}{'detail'}{'gnomad_exome'}{$1}{'freq'} = $arr[$h_header{$af_header}];
			}
			elsif( $af_header =~ /^gnomadGenome_(\S+)$/ ){
				$h_record{'popfreq'}{'detail'}{'gnomad_genome'}{$1}{'freq'} = $arr[$h_header{$af_header}];
			}
		}
		#大家都有的人群频率库
	}
	
	$h_record{'popfreq'}{'used_freq'}{"freq"} = -1;
	$h_record{'popfreq'}{'used_freq'}{"source"} = "not_detected";
	$h_record{'popfreq'}{'used_freq'}{"rule"} = "max";

	foreach my $t_db( keys(%{ $h_record{'popfreq'}{'detail'} }) ){
		foreach my $t_pop( keys(%{ $h_record{'popfreq'}{'detail'}{$t_db} }) ){
			#change . to -1 for population freq
			if( !defined($h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'}) || $h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'} eq '.' ){
				$h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'} = -1;
			}
			#get the max freq
			if( $h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'} > $h_record{'popfreq'}{'used_freq'}{"freq"} ){
				$h_record{'popfreq'}{'used_freq'}{"freq"} = $h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'freq'};
				$h_record{'popfreq'}{'used_freq'}{"source"} = "$t_db:$t_pop";
			}
			#reserve for hom/het count
			$h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'hom'} = -1;
			$h_record{'popfreq'}{'detail'}{$t_db}{$t_pop}{'het'} = -1;
		}
	}

	#人群频率标签
	$h_record{'tags'}{'freq_1'} = 1 if( $h_record{'popfreq'}{'used_freq'}{"freq"} <=0.01 );
	$h_record{'tags'}{'freq_5'} = 1 if( $h_record{'popfreq'}{'used_freq'}{"freq"} <=0.05 );
	
	#基于人群频率的证据项判断
	if( $h_record{'popfreq'}{'used_freq'}{"freq"} > 0.05 ){
		$h_record{'auto_acmg'}{'BA1'} = 'A';
	}
	elsif( $h_record{'popfreq'}{'used_freq'}{"freq"} > 0.015 ){
		$h_record{'auto_acmg'}{'BS1'} = 'S';
	}
	elsif( $h_record{'popfreq'}{'used_freq'}{"freq"} <= 0.005 ){
		$h_record{'auto_acmg'}{'PM2'} = 'M';
	}

	#预测
	#$h_record{'predicted'}{'note'} = encode("utf-8",decode("GB2312","PP3:多种统计方法预测出该变异会对基因或基因产物造成有害的影响，包括保守性预测、进化预测、剪接位点影响等\nBP4:多种统计方法预测出该变异会对基因或基因产物无影响，包括保守性预测、进化预测、剪接位点影响等。\nBP7:同义变异且预测不影响剪接。"));
	$h_record{'predicted'}{'note'} = decode_utf8("PP3:多种统计方法预测出该变异会对基因或基因产物造成有害的影响，包括保守性预测、进化预测、剪接位点影响等\nBP4:多种统计方法预测出该变异会对基因或基因产物无影响，包括保守性预测、进化预测、剪接位点影响等。\nBP7:同义变异且预测不影响剪接。");
	#my @arr_predict = ('REVEL','VEST4','SIFT','SIFT4G','Polyphen2_HDIV','Polyphen2_HVAR','LRT','MutationTaster','FATHMM','PROVEAN','MetaSVM','MetaLR','MetaRNN','M-CAP','PrimateAI','DEOGEN2');
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
	$h_predict{'综合评估'}{'MitoTIP'}{'note'} = '线粒体预测工具';
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
print "end at:\t".`date`;

#print join("\n",(sort {$a cmp $b} keys(%h_temp)))."\n";
#my $codec = BSON->new;
#my $doc = {'abc'=>"PP3:多种统计方法预测出该变异会对基因或基因产物造成有害的影响，包括保守性预测、进化预测、剪接位点影响等\nBP4:多种统计方法预测出该变异会对基因或基因产物无影响，包括保守性预测、进化预测、剪接位点影>响等。\nBP7:同义变异且预测不影响剪接。", };

#$collection->insert_one( $codec->encode_one($doc) );
#print `date`;
#my @data = $collection->find({"task_id"=>'A10010101'})->all();
#print `date`;
#$h_in{'all_muts'} = \@muts;
#$collection->insert_one( \%h_in);
