#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use MongoDB;
use Encode;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(rmtree);
use Cwd qw(abs_path);


my $db_col_name = $ARGV[0];

my %h_ch_disease = ();
my %h_hpo2gene = ();
my $line;



my %h_modifier = ();
open MODI,"$Bin/hpo_modifier.txt" or die $!;
#类型en 类型    HPO     英文    中文
#death_age       死亡年龄        HP:0001522      Death in infancy        婴儿期死亡
while( $line = <MODI> ){
	next if($line =~ /^#/);
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_modifier{$arr[2]}{'class'} = $arr[0];
	$h_modifier{$arr[2]}{'full_en'} = $arr[3];
	$h_modifier{$arr[2]}{'full_ch'} = decode_utf8($arr[4]);
}
close MODI;


my $chpo_file = $ARGV[1];
open IN,$chpo_file or die $!;
while( $line = <IN> ){
	chomp $line;
	if( $line =~ /INSERT INTO `dbs_disease`/ ){
		$line =~ s/NULL,/'-',/g;
		my @arr = split(/', '/,$line);
		$arr[2] = "OMIM:".$arr[2];
		$h_ch_disease{$arr[2]}{'ch_name'} = decode_utf8($arr[5]) if($arr[5] =~ /\S/ && $arr[5] ne "-");
		$h_ch_disease{$arr[2]}{'ch_descr'} = decode_utf8($arr[13]) if($arr[13] =~ /\S/ && $arr[13] ne "-");
		$h_ch_disease{$arr[2]}{'prevalence'} = decode_utf8($arr[10]) if($arr[10] =~ /\S/ && $arr[10] ne "-");
	}
}
close IN;

my $genes_to_phenotype = $ARGV[2];
open IN,$genes_to_phenotype or die $!;
while( $line = <IN> ){
	chomp $line;
	next if( $line =~ /^#/ );
	my @arr = split(/\t/,$line);
	$h_hpo2gene{$arr[8]}{"gene"}{$arr[1]}++ if( defined($arr[8]) && $arr[8] =~ /\S/ );
}
close IN;


while( $line = <STDIN> ){
	chomp $line;
	next if( $line =~ /^#/ );
	my @arr = split(/\t/,$line);
	my $relate_type = $arr[10];
	if( $arr[2] eq "NOT" ){
		$relate_type = "not_phenotype";
	}
	elsif( $arr[10] eq "P" ){
		$relate_type = "with_phenotype";
	}
	elsif( $arr[10] eq "I" ){
		$relate_type = "inher_mode";
	}
	else{
		if( defined( $h_modifier{$arr[3]} ) ){
			$relate_type = $h_modifier{$arr[3]}{'class'};
		}
	}
	$h_hpo2gene{$arr[0]}{"en_name"} = $arr[1] if( $arr[1] =~ /\S/ );
	$h_hpo2gene{$arr[0]}{$relate_type}{$arr[3]}{'flag'} = 1;
	$h_hpo2gene{$arr[0]}{$relate_type}{$arr[3]}{'refs'}{$arr[4]}++ if($arr[4] =~ /\S/);
	if( $relate_type eq "with_phenotype" || $relate_type eq "not_phenotype" ){
		if($arr[6] =~ /\S/){
			$h_hpo2gene{$arr[0]}{$relate_type}{$arr[3]}{'onset'}{$arr[6]}{$arr[4]}++;
		}
		if($arr[7] =~ /\S/){
			$h_hpo2gene{$arr[0]}{$relate_type}{$arr[3]}{'freq'}{$arr[7]}{$arr[4]}++;
		}
		if($arr[8] =~ /\S/){
			$h_hpo2gene{$arr[0]}{$relate_type}{$arr[3]}{'sex'}{$arr[8]}{$arr[4]}++;
		}
	}
}


my $client = MongoDB::MongoClient->new(
   host => '10.10.9.35:31112',
   username => 'lintop',
   password => 'lintop.hx321.mongo'
);
$client->reconnect;


my $collection = $client->ns( $db_col_name );
my %h_inher_ch = ();
$h_inher_ch{'HP:0000006'}{'full_en'} = 'Autosomal dominant inheritance';
$h_inher_ch{'HP:0000007'}{'full_en'} = 'Autosomal recessive inheritance';
$h_inher_ch{'HP:0001417'}{'full_en'} = 'X-linked inheritance';
$h_inher_ch{'HP:0001419'}{'full_en'} = 'X-linked recessive inheritance';
$h_inher_ch{'HP:0001423'}{'full_en'} = 'X-linked dominant inheritance';
$h_inher_ch{'HP:0001450'}{'full_en'} = 'Y-linked inheritance';

$h_inher_ch{'HP:0000006'}{'short_en'} = 'AD';
$h_inher_ch{'HP:0000007'}{'short_en'} = 'AR';
$h_inher_ch{'HP:0001417'}{'short_en'} = 'XL';
$h_inher_ch{'HP:0001419'}{'short_en'} = 'XLR';
$h_inher_ch{'HP:0001423'}{'short_en'} = 'XLD';
$h_inher_ch{'HP:0001450'}{'short_en'} = 'YL';

$h_inher_ch{'HP:0000006'}{'full_ch'} = decode_utf8('常染色体显性遗传');
$h_inher_ch{'HP:0000007'}{'full_ch'} = decode_utf8('常染色体隐性遗传');
$h_inher_ch{'HP:0001417'}{'full_ch'} = decode_utf8('X连锁遗传');
$h_inher_ch{'HP:0001419'}{'full_ch'} = decode_utf8('X连锁隐性遗传');
$h_inher_ch{'HP:0001423'}{'full_ch'} = decode_utf8('X连锁显性遗传');
$h_inher_ch{'HP:0001450'}{'full_ch'} = decode_utf8('Y连锁遗传');


$h_inher_ch{'HP:0000006'}{'short_ch'} = decode_utf8('常显');
$h_inher_ch{'HP:0000007'}{'short_ch'} = decode_utf8('常隐');
$h_inher_ch{'HP:0001417'}{'short_ch'} = decode_utf8('X连锁');
$h_inher_ch{'HP:0001419'}{'short_ch'} = decode_utf8('X隐');
$h_inher_ch{'HP:0001423'}{'short_ch'} = decode_utf8('X显');
$h_inher_ch{'HP:0001450'}{'short_ch'} = decode_utf8('Y连锁');
open MODE,"$Bin/hpo_inher_mode.txt" or die $!;
#hpo编号        英文缩写        中文缩写        中文名  英文名
#HP:0000006      AD      常显    常染色体显性遗传        Autosomal dominant inheritance
while( $line = <MODE> ){
	next if($line =~ /^#/);
	chomp $line;
	my @arr = split(/\t/,$line);
	$h_inher_ch{$arr[0]}{'short_en'} = $arr[1];
	$h_inher_ch{$arr[0]}{'short_ch'} = decode_utf8($arr[2]);
	$h_inher_ch{$arr[0]}{'full_ch'} = decode_utf8($arr[3]);
	$h_inher_ch{$arr[0]}{'full_en'} = $arr[4];
}
close MODE;


foreach my $dis( sort {$a cmp $b} keys(%h_hpo2gene) ){

	my %h_record= ();
	if( $dis =~ /OMIM/ ){
		$h_record{'type'} = "OMIM";
	}
	elsif( $dis =~ /ORPHA/ ){
		$h_record{'type'} = "ORPHA";
	}
	elsif( $dis =~ /DECIPHER/ ){
		$h_record{'type'} = "DECIPHER";
	}
	else{
		$h_record{'type'} = "O";
	}
	$h_record{'dis_id'} = "$dis";
	#$h_record{'dis_id'} =~ s/:/_/g;
	$h_record{'en_name'} = $h_hpo2gene{$dis}{'en_name'} if( defined( $h_hpo2gene{$dis}{'en_name'} ) );
	if( defined( $h_ch_disease{$dis}{'ch_name'} ) ){
		$h_record{'ch_name'} = $h_ch_disease{$dis}{'ch_name'};
	}
	if( defined( $h_ch_disease{$dis}{'prevalence'} ) ){
		my @arr_preva = ();
		my %h_preva = ();
		$h_preva{'ch_name'} = $h_ch_disease{$dis}{'prevalence'};
		$h_preva{'ref'} = "sakura";
		push@arr_preva,\%h_preva;
		$h_record{'prevalence'} = \@arr_preva;
		#$h_record{'ch_name'} = $h_ch_disease{$dis}{'prevalence'};
	}

	if( defined( $h_ch_disease{$dis}{'ch_descr'} ) ){
		$h_record{'ch_descr'} = $h_ch_disease{$dis}{'ch_descr'};
	}

	foreach my $gene(sort {$a cmp $b} keys(%{$h_hpo2gene{$dis}{"gene"}})){
		$h_record{'gene'}{$gene}{'name'} = $gene;
	}

	foreach my $relate_type("with_phenotype", "not_phenotype", "inher_mode", "onset_age", "death_age" ){
		foreach my $hpo(sort {$a cmp $b} keys(%{$h_hpo2gene{$dis}{$relate_type}})){
			my $p_id = $hpo;
			$p_id =~ s/:/_/;
			$h_record{$relate_type}{$p_id}{'hpo_id'} = $hpo;
			if( $relate_type eq "inher_mode" ){
				foreach my $name_type('full_ch','short_ch','short_en','full_en'){
					$h_record{$relate_type}{$p_id}{$name_type} = $hpo;
					if( defined($h_inher_ch{$hpo}{$name_type}) ){
						$h_record{$relate_type}{$p_id}{$name_type} = $h_inher_ch{$hpo}{$name_type};
					}
				}
			}
			elsif( $relate_type eq "onset_age" ||  $relate_type eq "death_age" ){
				foreach my $name_type('full_ch','full_en'){
					$h_record{$relate_type}{$p_id}{$name_type} = $hpo;
					if( defined($h_modifier{$hpo}{$name_type}) ){
						$h_record{$relate_type}{$p_id}{$name_type} = $h_modifier{$hpo}{$name_type};
						#print STDERR $h_record{$relate_type}{$p_id}{$name_type}."\n";
					}
				}
			}


			foreach my $t1('freq','sex','onset'){
				my @arr_t1 = ();
				foreach my $t2(sort {$a cmp $b} keys(%{$h_hpo2gene{$dis}{$relate_type}{$hpo}{$t1}})){
					my %h_t2 = ();
					$h_t2{'name'} = $t2;
					my @arr_temp_ref = ();
					foreach my $t_ref(sort {$a cmp $b} keys(%{$h_hpo2gene{$dis}{$relate_type}{$hpo}{$t1}{$t2}})){
						my %h_temp_ref = ();
						$h_temp_ref{'ref_id'} = $t_ref;
						push@arr_temp_ref,\%h_temp_ref;
					}
					$h_t2{'refs'} = \@arr_temp_ref;
					push @arr_t1,\%h_t2;
				}
				$h_record{$relate_type}{$p_id}{$t1} = \@arr_t1 if( scalar@arr_t1 > 0 );
			}
			#$h_record{$relate_type}{$p_id}{'freq'} = $h_hpo2gene{$dis}{$relate_type}{$hpo}{'freq'} if( defined($h_hpo2gene{$dis}{$relate_type}{$hpo}{'freq'}) );
			#$h_record{$relate_type}{$p_id}{'sex'}  = $h_hpo2gene{$dis}{$relate_type}{$hpo}{'sex'}  if( defined($h_hpo2gene{$dis}{$relate_type}{$hpo}{'sex'}) );
			if( defined($h_hpo2gene{$dis}{$relate_type}{$hpo}{'refs'}) ){
				my @arr_temp_ref = ();
				foreach my $t_ref(sort {$a cmp $b} keys(%{$h_hpo2gene{$dis}{$relate_type}{$hpo}{'refs'}})){
					my %h_temp_ref = ();
					$h_temp_ref{'ref_id'} = $t_ref;
					push@arr_temp_ref,\%h_temp_ref;
				}
				$h_record{$relate_type}{$p_id}{'refs'} = \@arr_temp_ref if(scalar@arr_temp_ref > 0);
			}
		}
	}

	$collection->insert_one( \%h_record );
}

$client->disconnect;
