#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.0";
my $DATE = "2022-03-08";
###################################################################

my($report, $mpa, $db, $name, $validate, $highconf, $type, $outdir, $version);
GetOptions(
	"report:s" => \$report,
	"mpa:s" => \$mpa,
	"db:s" => \$db,
	"name:s" => \$name,
	"validate:s" => \$validate,
	"highconf:s" => \$highconf,
	"type:s" => \$type,
	"outdir:s" => \$outdir,
	"version|v" => \$version,
);

if($version){
    print basename $0." $VERSION\n";
    exit(0);
}

&help unless($report && $mpa && $db && $validate && $highconf);
$name ||= "demo";
$type ||= "Others";
$outdir ||= "./";
$outdir = abs_path $outdir;
`mkdir -p $outdir` unless(-d $outdir);

my %column = (
	"RT" => 15,
	"Plasma" => 16,
	"CSF" => 17,
	"Others" => 18,
);

open IN,"$db/all.latin2chinese.INFO.txt" or die $!;
my %anno;
<IN>;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$anno{$t[2]}{info} = "$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]";
	$anno{$t[2]}{mark} = "$t[8]\t$t[9]";
}
close IN;

### format of the "$db/GroupList.info" is as follow:                                                                                                 ###
### GroupTaxid\tGroupName\tSpeciesTaxid\tSpeciesName\tGroupChineseName\tGenusTaxid\tGenus\tSuperkingdom\tFocus\Info\tReference\t[G+, G-]\tType_Judge ###
my (%group2show, %group2infos);
open IN,"$db/GroupList.info" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	#$anno{$t[1]}{info} = "$t[4]\t$t[7]\t$t[10]\t-\t$t[11]";
	#$anno{$t[1]}{mark} = "$t[12]\t$t[13]";
	$group2show{$t[1]} = 1 if($t[9] eq "*");
	$group2infos{$t[1]}{GroupTaxid} = $t[0];
	$group2infos{$t[1]}{GroupCN} = $t[4];
	$group2infos{$t[1]}{GenusCN} = $t[7];
	$group2infos{$t[1]}{Info} = $t[10];
	$group2infos{$t[1]}{Reference} = $t[11];
	$group2infos{$t[1]}{Genus} = $t[6];
	$group2infos{$t[1]}{Class} = $t[8];
	$group2infos{$t[1]}{TypeJudge} = $t[13];
}
close IN;

### format of the "$db/SpeciesGroup.list" is as follow:     ###
### GroupTaxid\tGroupName\tSeqTaxid\tSeqName                ###
open IN,"$db/SpeciesGroup.list" or die $!;
<IN>;
my(%name2group, %count);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$name2group{$t[3]} = $t[1];
	$count{$t[1]} = 1;
}
close IN;

my @groups = keys %count;

my @dbFile = ("$db/all.Taxonomy.txt", "$db/all.Taxonomy.other.txt");
my(%spelength, %strlength, %patho, %class, %trans, %speTaxid, %strTaxid, %name2genus, %taxid2spe, %target, %groupTran);
foreach my $i(@dbFile){
	open IN,"$i" or die $!;
	<IN>;
	while(<IN>){
		chomp;
		my @t = split /\t/,$_;
		$class{$t[8]} = $t[2];      # spe2class
		$class{$t[9]} = $t[2];
		$speTaxid{$t[8]} = $t[1];   # spe2spe_tid
		$strTaxid{$t[9]} = $t[0];   # strain2seq_tid
		$strlength{$t[9]} = $t[10]; # spe2strain_len
		$spelength{$t[8]} = $t[11]; # spe2spe_len
		$trans{$t[8]} = $t[14];   # spe2CN
		$trans{$t[9]} = $t[14];
		$name2genus{$t[8]} = $t[7]; # spe2genus
		$name2genus{$t[9]} = $t[7];
		if($t[2] eq 'Viruses'){
			$taxid2spe{$t[0]} = $t[9]; #seq_tid2spe
			$taxid2spe{$t[1]} = $t[8]; #spe_tid2spe
		}else{
			$taxid2spe{$t[0]} = $t[8];
			$taxid2spe{$t[1]} = $t[8];
		}
		$target{$t[8]} = $t[$column{$type}]; # spe2target
		$target{$t[9]} = $t[$column{$type}];
		if(exists $name2group{$t[8]}){
			$groupTran{$name2group{$t[8]}}{k} = $t[2];   # KingdomName
			$groupTran{$name2group{$t[8]}}{g} = $t[7];   # GenusName
		}
		if($t[12] eq '*'){
			$patho{$t[8]} = 1;
			$patho{$t[9]} = 1;
		}
	}
	close IN;
}

open IN,"$validate" or die $!;
my (%nodelete, %nodelete_groups, %check);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if($t[2] eq 'true'){
		$nodelete{$t[1]} = 1;
		$check{$t[1]} = "True";
		$nodelete_groups{$name2group{$t[1]}} = 1;
	}
}
close IN;

#####       Rewrote by Sujiawei at 2022/2/11      #####
open IN,"$validate" or die $!;
my (%delete, %deleted_groups);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if(!exists $nodelete{$t[1]}){
		$delete{$t[1]} = 1;
		if(!exists $nodelete_groups{$name2group{$t[1]}}){
			$deleted_groups{$name2group{$t[1]}} = 1;
		}
	}
}
close IN;

# one or more speicies under the group is false negative
# minus its/their reads count from the group reads count 
my (%group2minus, $tmp_group_name);
open IN,"$report" or die $!;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$t[5] =~ s/^\s+//g;
	if(exists $delete{$t[5]}){
		if(exists $name2group{$t[5]}){
			$tmp_group_name = $name2group{$t[5]};
			if(!exists $group2minus{$tmp_group_name}){
				$group2minus{$tmp_group_name} = $t[1];
			}else{
				$group2minus{$tmp_group_name} += $t[1];
			}
		}
	}
}

foreach my $spe(keys %delete){
	delete $name2group{$spe};
}
#######################################################

open IN,"$highconf" or die $!;
my %score;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$score{$taxid2spe{$t[0]}} = $t[1];
}
close IN;

open IN,"$mpa" or die $!;
my %stat = (
	"Viruses" => 0,            "Bacteria" => 0,
	"Archaea" => 0,            "Fungi" => 0,
	"Metazoa" => 0,            "Eukaryota" => 0,
	"Human" => 0,              "Unclassified" => 0,
	"Classified" => 0,         "Total_Reads" => 0,
	"Unclassified_Rate" => 0,  "Classified_Rate" => 0,
	"Metazoa_Parasite" => 0,   "Protozoa" => 0,
);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if($t[0] eq 'd__Viruses'){
		$stat{Viruses} = $t[1];
	}
	if($t[0] eq 'd__Bacteria'){
		$stat{Bacteria} = $t[1];
	}
	if($t[0] eq 'd__Archaea'){
		$stat{Archaea} = $t[1];
	}
	if($t[0] eq 'd__Eukaryota|k__Fungi'){
		$stat{Fungi} = $t[1];
	}
	if($t[0] eq 'd__Eukaryota|k__Metazoa'){
		$stat{Metazoa} = $t[1];
	}
	if($t[0] eq 'd__Eukaryota'){
		$stat{Eukaryota} = $t[1];
	}
	if($t[0] =~ /\|g__Homo$/){
		$stat{Human} = $t[1];
	}
}
close IN;

open IN,"$report" or die $!;
open TAX,">$outdir/taxid.list" or die $!;
my(%SpeLv, %sum, $preClass, $preSpe, %types, %reads, %genusMax);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$t[5] =~ s/^\s+//g;
	if($t[5] eq 'unclassified'){
		$stat{Unclassified} = $t[1];
	}
	if($t[5] eq 'root'){
		$stat{Classified} = $t[1];
	}
	# name => Reads_count
	$reads{$t[5]} = $t[1];
#	if(exists $delete{$t[5]} && exists $name2group{$t[5]}){
#		delete $reads{$name2group{$t[5]}};
#	}
#	if(!exists $nodelete{$t[5]} && exists $name2group{$t[5]}){
#		delete $reads{$name2group{$t[5]}};
#	}
	if($t[3] eq 'S' && $t[4] != 9606 && !exists $delete{$t[5]}){
		$preClass = $class{$t[5]};
		$preSpe = $t[5];
		if($preClass ne 'Viruses'){
			$SpeLv{$t[5]} = $t[1];
		}else{
			$types{$t[5]}{common} = $t[2];
		}
		$sum{S} += $t[1];
		print TAX "\n$preClass\t$preSpe\t$t[4]" if(exists $patho{$preSpe});
	}
	if($t[3] =~ /^S\d/ && exists $patho{$preSpe} && !exists $delete{$preSpe}){
		print TAX ",$t[4]";
	}
	if($preClass eq 'Viruses' && $t[3] eq 'S1'){
		$t[5] =~ s/^\s+//g;
		$types{$preSpe}{uniqTotal} += $t[1];
		$types{$preSpe}{types}{$t[5]} = $t[1];
	}
	
	#### Modified by Sujiawei at 2022/1/25 #####
	if($preClass eq 'Viruses' && $t[3] eq 'S2'){
		$t[5] =~ s/^\s+//g;
		if(exists $delete{$t[5]}){
			$types{$preSpe}{uniqTotal} -= $t[1];
			if($types{$preSpe}{uniqTotal} eq 0){
				delete $types{$preSpe};
			}
		}
	}
    ############################################

	if(exists $name2genus{$t[5]} && $class{$t[5]} ne 'Viruses'){
		$genusMax{$name2genus{$t[5]}}{$class{$t[5]}} = $t[1] if(!exists $genusMax{$name2genus{$t[5]}}{$class{$t[5]}} || $genusMax{$name2genus{$t[5]}}{$class{$t[5]}} < $t[1]);
	}
}
close IN;
close TAX;

#####       Rewrote by Sujiawei at 2022/2/11      #####
my (%existGroup, %group2reads_cnt);
foreach my $i(@groups){
	if(exists $reads{$i}){
		$existGroup{$i} = $reads{$i};
		$group2reads_cnt{$i} = $reads{$i};
	}
}
#######################################################

my %strain2species;
foreach my $i(keys %types){
	my $number =  scalar (keys %{$types{$i}{types}});
	if($number == 0){
		$SpeLv{$i} = $types{$i}{common};
		if(exists $name2genus{$i}){
			$genusMax{$name2genus{$i}}{$class{$i}} = $SpeLv{$i} if(!exists $genusMax{$name2genus{$i}}{$class{$i}} || $genusMax{$name2genus{$i}}{$class{$i}} < $SpeLv{$i});
		}
	}else{
		foreach my $j(keys %{$types{$i}{types}}){
			$SpeLv{$j} = sprintf("%d", $types{$i}{types}{$j} / $types{$i}{uniqTotal} * $types{$i}{common} + $types{$i}{types}{$j});
			if(exists $name2genus{$j}){
				$genusMax{$name2genus{$j}}{$class{$j}} = $SpeLv{$j} if(!exists $genusMax{$name2genus{$j}}{$class{$j}} || $genusMax{$name2genus{$j}}{$class{$j}} < $SpeLv{$j});
			}
			$strain2species{$j} = $i;
		}
	}
}

my(%uniform, %RPM);
foreach my $sp(keys %SpeLv){
	my $value;
	if(exists $strain2species{$sp}){
		if(exists $strlength{$sp}){
			$value = $SpeLv{$sp} / $strlength{$sp};
		}else{
			$value = $SpeLv{$sp} / $spelength{$strain2species{$sp}};
		}
	}else{
		$value = $SpeLv{$sp} / $spelength{$sp} if($spelength{$sp});
	}
	if(!exists $delete{$sp}){
		$uniform{$sp} = $value;
		$sum{uniform} += $value;
		$RPM{$sp} = sprintf("%.2f", $SpeLv{$sp} / $sum{S} * 1000000);
	}
}

##### Added by Sujiawei at 2022/2/11 #####
foreach my $gp_name(keys %group2minus){
	if(exists $group2reads_cnt{$gp_name}){
		$group2reads_cnt{$gp_name} -= $group2minus{$gp_name};
	}
}
############################################

open OU1,">$outdir/Taxonomy_Summary.txt" or die $!;
print OU1 "Sample\tTotal_Reads\tUnclassified_Reads_Number\tUnclassified_Rate\tClassified_Reads_Number\tClassified_Rate\tViruses\tBacteria\tArchaea\tFungi\tProtozoa\tMetazoa_Parasite\tHuman\n";
$stat{Total_Reads} = $stat{Unclassified} + $stat{Classified};
$stat{Unclassified_Rate} = sprintf("%.2f", $stat{Unclassified}/$stat{Total_Reads}*100);
$stat{Classified_Rate} = sprintf("%.2f", $stat{Classified}/$stat{Total_Reads}*100);
$stat{Metazoa_Parasite} = $stat{Metazoa} - $stat{Human};
$stat{Protozoa} = $stat{Eukaryota} - $stat{Fungi} - $stat{Metazoa};
print OU1 "$name\t$stat{Total_Reads}\t$stat{Unclassified}\t$stat{Unclassified_Rate}\t$stat{Classified}\t$stat{Classified_Rate}";
print OU1 "\t$stat{Viruses}\t$stat{Bacteria}\t$stat{Archaea}\t$stat{Fungi}\t$stat{Protozoa}\t$stat{Metazoa_Parasite}\t$stat{Human}\n";
close OU1;

open OU2,">$outdir/Total_Detail.txt" or die $!;
open OU3,">$outdir/Pathogeny_Detail.txt" or die $!;
open OU4,">$outdir/bacteria.xls" or die $!;
open OU5,">$outdir/fungi.xls" or die $!;
open OU6,">$outdir/LY.xls" or die $!;
open OU7,">$outdir/protozoa.xls" or die $!;
open OU8,">$outdir/viral.xls" or die $!;
open IN,"$outdir/taxid.list" or die $!;
open TMP,">$outdir/taxid.list.tmp" or die $!;
print OU2 "Kingdom\tGenusName\tGenusReads\tGroupName\tGroupReads\tScientificName\tChineseName\tReads_Number\tGMRN\tTPM\tRH-RPM\tVerification\tFocus\n";
print OU3 "Kingdom\tGenusName\tGenusReads\tGroupName\tGroupReads\tScientificName\tChineseName\tReads_Number\tGMRN\tTPM\tRH-RPM\tVerification\tFocus\n";
print OU4 "tax_id\tspecies\treads_count\tcoverage\tspecies_cn\tgenus_cn\tnoun\tmedicine\tncbi\tblast_maxscore\tmark\ttype_judge\n";
print OU5 "tax_id\tspecies\treads_count\tcoverage\tspecies_cn\tgenus_cn\tnoun\tmedicine\tncbi\tblast_maxscore\tmark\ttype_judge\n";
print OU6 "tax_id\tspecies\treads_count\tcoverage\tspecies_cn\tgenus_cn\tnoun\tmedicine\tncbi\tblast_maxscore\tmark\ttype_judge\n";
print OU7 "tax_id\tspecies\treads_count\tcoverage\tspecies_cn\tgenus_cn\tnoun\tmedicine\tncbi\tblast_maxscore\tmark\ttype_judge\n";
print OU8 "tax_id\tspecies\treads_count\tcoverage\tspecies_cn\tgenus_cn\tnoun\tmedicine\tncbi\tblast_maxscore\tmark\ttype_judge\n";
my(%sort, %save);
foreach my $sp(sort {sortByNameAfterNum(\%RPM)} keys %RPM){
	my $UV = sprintf("%.2f", $uniform{$sp} / $sum{uniform} * 1000000);
	if(!exists $anno{$sp}){
		$anno{$sp}{info} = "-\t-\t-\t-\t-";
		$anno{$sp}{mark} = "-\t-";
	}
	if(exists $class{$sp}){
		print OU2 "$class{$sp}";
		if($name2genus{$sp} eq 'NA'){
			print OU2 "\t-\t0";
		}else{
			print OU2 "\t$name2genus{$sp}\t$reads{$name2genus{$sp}}";
		}
		if(exists $name2group{$sp}){
			print OU2 "\t$name2group{$sp}\t$group2reads_cnt{$name2group{$sp}}";
		}else{
			print OU2 "\t-\t0";
		}
		print OU2 "\t$sp\t$trans{$sp}\t$SpeLv{$sp}\t$genusMax{$name2genus{$sp}}{$class{$sp}}\t$UV\t$RPM{$sp}";
		if(exists $check{$sp}){
			print OU2 "\t$check{$sp}\t$target{$sp}\n";
		}else{
			if(exists $score{$sp}){
				print OU2 "\t$score{$sp}\t$target{$sp}\n";
			}else{
				print OU2 "\tFalse\t$target{$sp}\n";
			}
		}
		
		if(exists $existGroup{$name2group{$sp}} && !exists $group2show{$name2group{$sp}}){
			delete $existGroup{$name2group{$sp}};
		}
	}
	if(exists $patho{$sp} && exists $class{$sp}){
		if($class{$sp} ne 'Viruses'){
			unless($sp =~ / sp\.| genomosp\.| genosp\./){
				print OU3 "$class{$sp}";
				if($name2genus{$sp} eq 'NA'){
					print OU3 "\t$name2genus{$sp}\t0";
				}else{
					print OU3 "\t$name2genus{$sp}\t$reads{$name2genus{$sp}}";
				}
				if(exists $name2group{$sp}){
					print OU3 "\t$name2group{$sp}\t$group2reads_cnt{$name2group{$sp}}";
				}else{
					print OU3 "\t-\t0";
				}
				print OU3 "\t$sp\t$trans{$sp}\t$SpeLv{$sp}\t$genusMax{$name2genus{$sp}}{$class{$sp}}\t$UV\t$RPM{$sp}";
				if(exists $check{$sp}){
					print OU3 "\t$check{$sp}\t$target{$sp}\n";
				}else{
					if(exists $score{$sp}){
						print OU3 "\t$score{$sp}\t$target{$sp}\n";
					}else{
						print OU3 "\t-\t$target{$sp}\n";
					}
				}
			}
			if($class{$sp} eq 'Bacteria'){
				print OU4 "$speTaxid{$sp}\t$sp\t$SpeLv{$sp}\t-\t$anno{$sp}{info}";
				if(exists $check{$sp}){
					print OU4 "\t$check{$sp}\t$anno{$sp}{mark}\n";
				}else{
					if(exists $score{$sp}){
						print OU4 "\t$score{$sp}\t$anno{$sp}{mark}\n";
					}else{
						print OU4 "\t-\t$anno{$sp}{mark}\n";
					}
				}
			}elsif($class{$sp} eq 'Eukaryota:Fungi'){
				print OU5 "$speTaxid{$sp}\t$sp\t$SpeLv{$sp}\t-\t$anno{$sp}{info}";
				if(exists $check{$sp}){
					print OU5 "\t$check{$sp}\t$anno{$sp}{mark}\n";
				}else{
					if(exists $score{$sp}){
						print OU5 "\t$score{$sp}\t$anno{$sp}{mark}\n";
					}else{
						print OU5 "\t-\t$anno{$sp}{mark}\n";
					}
				}
			}elsif($class{$sp} eq 'Eukaryota:Protozoa' || $class{$sp} eq 'Eukaryota:Parasite'){
				print OU7 "$speTaxid{$sp}\t$sp\t$SpeLv{$sp}\t-\t$anno{$sp}{info}";
				if(exists $check{$sp}){
					print OU7 "\t$check{$sp}\t$anno{$sp}{mark}\n";
				}else{
					if(exists $score{$sp}){
						print OU7 "\t$score{$sp}\t$anno{$sp}{mark}\n";
					}else{
						print OU7 "\t-\t$anno{$sp}{mark}\n";
					}
				}
			}else{
				print OU6 "$speTaxid{$sp}\t$sp\t$SpeLv{$sp}\t-\t$anno{$sp}{info}";
				if(exists $check{$sp}){
					print OU6 "\t$check{$sp}\t$anno{$sp}{mark}\n";
				}else{
					if(exists $score{$sp}){
						print OU6 "\t$score{$sp}\t$anno{$sp}{mark}\n";
					}else{
						print OU6 "\t-\t$anno{$sp}{mark}\n";
					}
				}
			}
		}else{
			unless($sp =~ / sp\.| genomosp\.| genosp\./){
				print OU3 "$class{$sp}";
				if($name2genus{$sp} eq 'NA'){
					print OU3 "\t$name2genus{$sp}\t0";
				}else{
					print OU3 "\t$name2genus{$sp}\t$reads{$name2genus{$sp}}";
				}
				if(exists $name2group{$sp}){
					print OU3 "\t$name2group{$sp}\t$group2reads_cnt{$name2group{$sp}}";
				}else{
					print OU3 "\t-\t0";
				}
				print OU3 "\t$sp\t$trans{$sp}\t$SpeLv{$sp}\t$genusMax{$name2genus{$sp}}{$class{$sp}}\t$UV\t$RPM{$sp}";
				if(exists $check{$sp}){
					print OU3 "\t$check{$sp}\t$target{$sp}\n";
				}else{
					if(exists $score{$sp}){
						print OU3 "\t$score{$sp}\t$target{$sp}\n";
					}else{
						print OU3 "\tFalse\t$target{$sp}\n";
					}
				}
				# output to viral.xls
				if(exists $check{$sp}){
					print OU8 "$strTaxid{$sp}\t$sp\t$SpeLv{$sp}\t-\t$anno{$sp}{info}\t$check{$sp}\t$anno{$sp}{mark}\n";
				}else{
					if(exists $score{$sp}){
						print OU8 "$strTaxid{$sp}\t$sp\t$SpeLv{$sp}\t-\t$anno{$sp}{info}\t$score{$sp}\t$anno{$sp}{mark}\n";
					}else{
						print OU8 "$strTaxid{$sp}\t$sp\t$SpeLv{$sp}\t-\t$anno{$sp}{info}\tFalse\t$anno{$sp}{mark}\n"
					}
				}
			}
		}
		if(exists $strain2species{$sp}){
			$sort{$class{$sp}}{$strain2species{$sp}} = 1;
			$sort{$class{$sp}}{$strain2species{$sp}} = scalar %{$sort{$class{$sp}}};
		}else{
			$sort{$class{$sp}}{$sp} = 1;
			$sort{$class{$sp}}{$sp} = scalar %{$sort{$class{$sp}}};
		}
	}
}

##### Rewrote by Sujiawei at 2022/2/11 #####
foreach my $del_group(keys %deleted_groups){
	delete $existGroup{$del_group} if(exists $existGroup{$del_group})
}

foreach my $gp_name(keys %group2minus){
	if(exists $existGroup{$gp_name}){
		$existGroup{$gp_name} -= $group2minus{$gp_name};
	}
}

my $row;
foreach my $i(keys %existGroup){
	next if(!exists $group2show{$i});
	$row = "$groupTran{$i}{k}\t$groupTran{$i}{g}";
	if(exists $reads{$groupTran{$i}{g}}){
		$row .= "\t$reads{$groupTran{$i}{g}}";
	}else{
		$row .= "\t0";
	}
	$row .= "\t$i\t$existGroup{$i}\t$i\t$group2infos{$i}{GroupCN}\t$existGroup{$i}\t$genusMax{$group2infos{$i}{Genus}}{$group2infos{$i}{Class}}\t-\t-\t-\t*\n";
}

print OU2 "$row";
print OU3 "$row";

# add data of $existGroup to *xls
foreach my $g(keys %existGroup){
	if(!exists $group2show{$g}){
		next;
	}
	if($group2infos{$g}{Class} eq 'Viruses'){
		print OU8 "$group2infos{$g}{GroupTaxid}\t$g\t$existGroup{$g}\t-\t$group2infos{$g}{GroupCN}\t$group2infos{$g}{GenusCN}\t$group2infos{$g}{Info}\t-\t$group2infos{$g}{Reference}\t-\t-\t$group2infos{$g}{TypeJudge}\n";
	}else{
		if($group2infos{$g}{Class} eq 'Bacteria'){
			print OU4 "$group2infos{$g}{GroupTaxid}\t$g\t$existGroup{$g}\t-\t$group2infos{$g}{GroupCN}\t$group2infos{$g}{GenusCN}\t$group2infos{$g}{Info}\t-\t$group2infos{$g}{Reference}\t-\t-\t$group2infos{$g}{TypeJudge}\n";
		}elsif($group2infos{$g}{Class} eq 'Eukaryota:Fungi'){
			print OU5 "$group2infos{$g}{GroupTaxid}\t$g\t$existGroup{$g}\t-\t$group2infos{$g}{GroupCN}\t$group2infos{$g}{GenusCN}\t$group2infos{$g}{Info}\t-\t$group2infos{$g}{Reference}\t-\t-\t$group2infos{$g}{TypeJudge}\n";
		}elsif($group2infos{$g}{Class} eq 'Eukaryota:Protozoa' || $group2infos{$g}{Class} eq 'Eukaryota:Protozoa'){
			print OU7 "$group2infos{$g}{GroupTaxid}\t$g\t$existGroup{$g}\t-\t$group2infos{$g}{GroupCN}\t$group2infos{$g}{GenusCN}\t$group2infos{$g}{Info}\t-\t$group2infos{$g}{Reference}\t-\t-\t$group2infos{$g}{TypeJudge}\n";
		}else{
			print OU6 "$group2infos{$g}{GroupTaxid}\t$g\t$existGroup{$g}\t-\t$group2infos{$g}{GroupCN}\t$group2infos{$g}{GenusCN}\t$group2infos{$g}{Info}\t-\t$group2infos{$g}{Reference}\t-\t-\t$group2infos{$g}{TypeJudge}\n";
		}
	}
}

###########################################

close OU2; close OU3; close OU4; close OU5;
close OU6; close OU7; close OU8;

while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	print TMP "$_\t$sort{$t[0]}{$t[1]}\n";
}
close IN; close TMP;
`mv $outdir/taxid.list.tmp $outdir/taxid.list`;

sub sortByNameAfterNum{
	my $ref_hash = shift;
	if($ref_hash->{$a} != $ref_hash->{$b}){
		return $ref_hash->{$b} <=> $ref_hash->{$a};
	}else{
		return $a cmp $b;
	}
}

sub help{
print "
	Usage: perl $0

	--report             kraken's report file
	--mpa                kraken's mpa-report file
	--db                 path of the database to annotate
	--name               sample name
	--validate           Mark2Fa.pl's output
	--highconf           High confidence kmer score
	--type               RT|Plasma|CSF|[Others]
	--outdir             output path. [./]
	--version            print version information.

";
exit(0);
}
