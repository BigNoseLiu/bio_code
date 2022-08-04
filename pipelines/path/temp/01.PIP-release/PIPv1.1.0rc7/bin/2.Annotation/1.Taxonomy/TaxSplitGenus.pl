#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use FindBin qw($Bin);
use lib "$Bin/../../../lib/perl5";
use MyModule::GlobalVar qw($PIGZ $SRC $BWA $SAMTOOLS $BLASTDB $BLASTN);
###################################################################
my $AUTHOR = "qiuwancen";
my $EMAIL = "972538446\@qq.com";
my $VERSION = "V1.0";
my $DATE = "2022-04-01";
###################################################################

my($report, $inMark, $db, $ntV, $fq, $krakenOut, $outdir, $version);
GetOptions(
	"report:s" => \$report,
	"inMark:s" => \$inMark,
	"db:s" => \$db,
	"ntV:s" => \$ntV,
	"fq:s" => \$fq,
	"krakenout:s" => \$krakenOut,
	"outdir:s" => \$outdir,
	"version|v" => \$version,
);
if($version){
    print basename $0." $VERSION\n";
    exit(0);
}
&help unless($report && $inMark && $db && $fq && $krakenOut);
$ntV ||= "20210103";
$outdir ||= "./";
$outdir = abs_path $outdir;
`mkdir -p $outdir` unless(-d $outdir);

############################### Annotating by DB ###############################
my @dbFile = ("$db/all.Taxonomy.txt", "$db/all.Taxonomy.other.txt");
my(%anno, %seqid2gName, %seqid2speid, %id2kingdom);
foreach my $file(@dbFile){
	open IN,"$file" or die $!;
	while(<IN>){
		chomp;
		my @t = split /\t/,$_;
		$seqid2gName{$t[0]} = $t[7];
		$seqid2speid{$t[0]} = $t[1];
		$id2kingdom{$t[0]} = $t[2];
		$id2kingdom{$t[1]} = $t[2];
		if($t[2] eq 'Viruses'){
			$anno{$t[0]} = $t[9];
			$anno{$t[1]} = $t[8];
		}else{
			$anno{$t[0]} = $t[8];
		}
		# annotate genus taxid
		$anno{$t[13]} = $t[7]
	}
}

open IN,"$db/pathoSeqTaxID2SpeTaxID2Tax.list" or die $!;
#my(%anno, %seqid2gName, %seqid2speid, %id2kingdom);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	next if(exists $anno{$t[0]});
	$seqid2gName{$t[0]} = $t[7];
	$seqid2speid{$t[0]} = $t[1];
	$id2kingdom{$t[0]} = $t[2];
	$id2kingdom{$t[1]} = $t[2];
	if($t[2] eq 'Viruses'){
		$anno{$t[0]} = $t[9];
		$anno{$t[1]} = $t[8];
	}else{
		$anno{$t[0]} = $t[8];
	}
}
close IN;

open IN,"$db/SeqTaxid2GenusTaxid-K.patho.list" or die $!;
my(%taxid, %gID2gName);
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$taxid{$t[1]}{kingdom} = $t[0];
	$taxid{$t[1]}{genus} = $t[2];
	$id2kingdom{$t[1]} = $t[0] if(!exists $id2kingdom{$t[1]});
	$id2kingdom{$t[2]} = $t[0] if(!exists $id2kingdom{$t[2]});
	$gID2gName{$t[0]}{$t[2]} = $seqid2gName{$t[1]} if(exists $seqid2gName{$t[1]});
}
close IN;

open IN,"$db/SpeciesGroup.list" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$anno{$t[0]} = $t[1];
	$id2kingdom{$t[0]} = $id2kingdom{$t[2]};
}
close IN;

open IN, "$db/GroupList.info" or die $!;
<IN>;
while(<IN>){
    chomp;
    my @t = split /\t/,$_;
    next if(exists $taxid{$t[0]});
    $taxid{$t[0]}{kingdom} = $t[8];
    $taxid{$t[0]}{genus} = $t[5];
    $gID2gName{$t[8]}{$t[5]} = $t[6];
}

############################### Annotating by Analysis Result ###############################
open IN,"$inMark" or die $!;
my %hash;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$t[1] =~ s/\/.*//g;
	$hash{$t[1]} = $t[2];
}
close IN;

open RP,"$report" or die $!;
my(%genus2kingdom, %id2genus);
while(<RP>){
	chomp;
	my @t = split /\t/,$_;
	if(exists $taxid{$t[4]}){
	    if($taxid{$t[4]}{genus} eq 'NULL'){
	        $genus2kingdom{"NULL_$taxid{$t[4]}{kingdom}"} = $taxid{$t[4]}{kingdom};
	        $id2genus{$t[4]} = "NULL_$taxid{$t[4]}{kingdom}";
	        $gID2gName{$taxid{$t[4]}{kingdom}}{"NULL_$taxid{$t[4]}{kingdom}"} = "NULL_$taxid{$t[4]}{kingdom}";
	    }else{
	        $genus2kingdom{$taxid{$t[4]}{genus}} = $taxid{$t[4]}{kingdom};
	        $id2genus{$t[4]} = $taxid{$t[4]}{genus};
	    }
	}
}
close RP;
############################### Extract Result ###############################
my @FQs = split /,/, $fq;
my @Genus = keys %genus2kingdom;

open OU,">$outdir/validated.fna" or die $!;
open FAB,"| $PIGZ -c -p 16 >$outdir/bacteria.result.fq.gz" or die $!;
open FAF,"| $PIGZ -c -p 16 >$outdir/fungi.result.fq.gz" or die $!;
open FAV,"| $PIGZ -c -p 16 >$outdir/viral.result.fq.gz" or die $!;
open FAP,"| $PIGZ -c -p 16 >$outdir/protozoa.result.fq.gz" or die $!;
open FALY,"| $PIGZ -c -p 16 >$outdir/LY.result.fq.gz" or die $!;
open CHECK_LOG,"| $PIGZ -c -p 16 >$outdir/SequencesID_need_to_check.log" or die $!;
my %fh;
if(scalar @FQs == 1){
	if($FQs[0] =~ /\.fq\.gz$/ || $FQs[0] =~ /\.fastq\.gz$/){
		open FQ,"$PIGZ -dc -p 16 $FQs[0] |" or die $!;
	}elsif($FQs[0] =~ /\.fq$/ || $FQs[0] =~ /\.fastq$/){
		open FQ,"$FQs[0]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	my $name = basename $FQs[0];    $name =~ s/\.gz//g;
#	open SH,">$outdir/BWAcmd.list" or die $!;
	for(@Genus){
		my $gPath = "$outdir/Fasta/$genus2kingdom{$_}/$gID2gName{$genus2kingdom{$_}}{$_}";
		`mkdir -p $gPath` unless(-d $gPath);
		open $fh{$_},"| $PIGZ -c -p 16 >$gPath/$_.fna.gz" or die $!;
#		my $refname;
#		if($_ eq 'NULL'){
#			$refname = "NULL_$genus2kingdom{$_}";
#		}else{
#			$refname = $_;
#		}
#		print SH "$genus2kingdom{$_}\t$gID2gName{$genus2kingdom{$_}}{$_}\t$outdir\tFasta/$genus2kingdom{$_}/$gID2gName{$genus2kingdom{$_}}{$_}/$_.fna.gz\t$BWA mem -t 1 -Y -M $db/BWA/$refname/$refname.fna $gPath/$_.fna.gz |$SAMTOOLS sort - |$SAMTOOLS view -F4 -bS -o $gPath/$_.bam -\n";
	}
#	close SH;
	open KO,"$krakenOut" or die $!;
	my($fqID, $fqSeq, $fqQual);
	my $fqID2 = 0;
	while(<KO>){
		chomp;
		my($Rid, $Tid) = (split /\t/,$_)[1,2];
		$Rid =~ s/\/.*//g;
		if($id2kingdom{$Tid} != 'Viruses'){
			next unless(exists $id2genus{$Tid});
		}
		while($fqID2 ne $Rid){
			chomp($fqID = <FQ>);
			chomp($fqSeq = <FQ>);
			<FQ>;
			chomp($fqQual = <FQ>);
			$fqID =~ s/^\@//g;
			$fqID2 = $fqID;
			$fqID2 =~ s/\/.*$//g;   $fqID2 =~ s/ .*//g;
		}
		if(exists $anno{$Tid}){
		    $fh{$id2genus{$Tid}}->print(">$fqID $anno{$Tid}\n$fqSeq\n") if(exists $id2genus{$Tid});
			if($id2kingdom{$Tid} eq 'Bacteria'){
				print FAB ">$fqID $anno{$Tid}\n$fqSeq\n";
			}elsif($id2kingdom{$Tid} eq 'Eukaryota:Fungi'){
				print FAF ">$fqID $anno{$Tid}\n$fqSeq\n";
			}elsif($id2kingdom{$Tid} eq 'Viruses'){
				print FAV ">$fqID $anno{$Tid}\n$fqSeq\n";
			}elsif($id2kingdom{$Tid} eq 'Eukaryota:Parasite' || $id2kingdom{$Tid} eq 'Eukaryota:Protozoa'){
				print FAP ">$fqID $anno{$Tid}\n$fqSeq\n";
			}elsif($id2kingdom{$Tid} eq 'SpecialPathogens'){
				print FALY ">$fqID $anno{$Tid}\n$fqSeq\n";
			}else{
				print CHECK_LOG "check: $fqID\n";
			}
		}else{
			$fh{$id2genus{$Tid}}->print(">$fqID\n$fqSeq\n") if(exists $id2genus{$Tid});
			if($id2kingdom{$Tid} eq 'Bacteria'){
				print FAB ">$fqID\n$fqSeq\n";
			}elsif($id2kingdom{$Tid} eq 'Eukaryota:Fungi'){
				print FAF ">$fqID\n$fqSeq\n";
			}elsif($id2kingdom{$Tid} eq 'Viruses'){
				print FAV ">$fqID\n$fqSeq\n";
			}elsif($id2kingdom{$Tid} eq 'Eukaryota:Parasite' || $id2kingdom{$Tid} eq 'Eukaryota:Protozoa'){
				print FAP ">$fqID\n$fqSeq\n";
			}elsif($id2kingdom{$Tid} eq 'SpecialPathogens'){
				print FALY ">$fqID\n$fqSeq\n";
			}else{
				print CHECK_LOG "check: $fqID\n";
			}
		}
		print OU ">$fqID|kraken:taxid|$Tid\n$fqSeq\n" if(exists $hash{$fqID2});
	}
	close KO;
	for(keys %genus2kingdom){
		close $fh{$_};
	}
}else{
	if($FQs[0] =~ /\.fq\.gz$/ || $FQs[0] =~ /\.fastq\.gz$/){
		open FQ1,"$PIGZ -dc -p 16 $FQs[0] |" or die $!;
	}elsif($FQs[0] =~ /\.fq$/ || $FQs[0] =~ /\.fastq$/){
		open FQ1,"$FQs[0]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	if($FQs[1] =~ /\.fq\.gz$/ || $FQs[1] =~ /\.fastq\.gz$/){
		open FQ2,"$PIGZ -dc -p 16 $FQs[1] |" or die $!;
	}elsif($FQs[1] =~ /\.fq$/ || $FQs[1] =~ /\.fastq$/){
		open FQ2,"$FQs[1]" or die $!;
	}else{
		die "Please confirm the format of input reads file!\n";
	}
	my $name1 = basename $FQs[0];   $name1 =~ s/\.gz//g;
	my $name2 = basename $FQs[1];   $name2 =~ s/\.gz//g;
	open SH,">$outdir/BWAcmd.list" or die $!;
	for(@Genus){
		my $gPath = "$outdir/Fasta/$genus2kingdom{$_}/$gID2gName{$genus2kingdom{$_}}{$_}";
		`mkdir -p $gPath` unless(-d $gPath);
		open $fh{$_}{1},"| $PIGZ -c -p 16 >$gPath/$_.1.fna.gz" or die $!;
		open $fh{$_}{2},"| $PIGZ -c -p 16 >$gPath/$_.2.fna.gz" or die $!;
		my $refname;
		if($_ eq 'NULL'){
			$refname = "NULL_$genus2kingdom{$_}";
		}else{
			$refname = $_;
		}
		print SH "$genus2kingdom{$_}\t$gID2gName{$genus2kingdom{$_}}{$_}\t$outdir\tFasta/$genus2kingdom{$_}/$gID2gName{$genus2kingdom{$_}}{$_}/$_.1.fna.gz,Fasta/$genus2kingdom{$_}/$gID2gName{$genus2kingdom{$_}}{$_}/$_.2.fna.gz\t$BWA mem -t 1 -Y -M $db/BWA/$refname/$refname.fna $gPath/$_.1.fna.gz $gPath/$_.2.fna.gz |$SAMTOOLS sort - |$SAMTOOLS view -F4 -bS -o $gPath/$_.bam -\n";
	}
	close SH;
	open KO,"$krakenOut" or die $!;
	my($fqID1, $fqID2, $fqID2_2, $fqSeq1, $fqSeq2, $fqQual1, $fqQual2);
	my $fqID1_2 = 0;
	while(<KO>){
		chomp;
		my($Rid, $Tid) = (split /\t/,$_)[1,2];
		if($id2kingdom{$Tid} != 'Viruses'){
			next unless(exists $id2genus{$Tid});
		}
		while($fqID1_2 ne $Rid){
			chomp($fqID1 = <FQ1>);    chomp($fqID2 = <FQ2>);
			chomp($fqSeq1 = <FQ1>);   chomp($fqSeq2 = <FQ2>);
			<FQ1>;                    <FQ2>;
			chomp($fqQual1 = <FQ1>);  chomp($fqQual2 = <FQ2>);
			$fqID1 =~ s/^\@//g;
			$fqID1_2 = $fqID1;  $fqID1_2 =~ s/\/.*$//g;   $fqID1_2 =~ s/ .*//g;
		}
		$fh{$id2genus{$Tid}}{1}->print(">$fqID1\n$fqSeq1\n") if(exists $id2genus{$Tid});
		$fh{$id2genus{$Tid}}{2}->print(">$fqID2\n$fqSeq2\n") if(exists $id2genus{$Tid});
		print OU ">$fqID1|kraken:taxid|$Tid\n$fqSeq1\n>$fqID2|kraken:taxid|$Tid\n$fqSeq2\n" if(exists $hash{$fqID1_2});
	}
	close KO;
	for(keys %genus2kingdom){
		close $fh{$_}{1};
		close $fh{$_}{2};
	}
}
close OU; close FAB; close FAF;
close FAV; close FAP; close FALY;

`$BLASTN -query $outdir/validated.fna -db $BLASTDB/$ntV/nt -num_threads 20 -out $outdir/validated.m8 -outfmt \"6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp\" -evalue 0.05 -word_size 28 -reward 1 -penalty -2 -perc_identity 95 -qcov_hsp_perc 90`; #-num_alignments 100000

open IN,"$outdir/validated.m8" or die $!;
my %score;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if(!exists $score{$t[0]} || $score{$t[0]}{"identity+coverage"} < $t[2]+$t[14]){
		$score{$t[0]}{"identity+coverage"} = $t[2] + $t[14];
	}
}
close IN;

open IN,"$outdir/validated.m8" or die $!;
my %mark;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	my $queryTaxID = (split /\|/,$t[0])[-1];
	my $targetTaxID;
	if($t[1] =~ /:/){
		$targetTaxID = (split /:/,$t[1])[0];
	}else{
		$targetTaxID = (split /\|/,$t[1])[0];
	}
	next unless($targetTaxID);
	if($t[2] + $t[14] == $score{$t[0]}{"identity+coverage"}){
		if(exists $seqid2speid{$queryTaxID} && exists $seqid2speid{$targetTaxID} && $seqid2speid{$queryTaxID} == $seqid2speid{$targetTaxID}){
			$mark{$t[0]} = 1;
		}
	}
}
close IN;

open VO,">$outdir/validated.out" or die $!;
foreach my $i(keys %score){
	my $queryTaxID = (split /\|/,$i)[-1];
	if(exists $mark{$i}){
		print VO "$i\t$anno{$queryTaxID}\ttrue\n" if(exists $anno{$queryTaxID});
	}else{
		print VO "$i\t$anno{$queryTaxID}\tfalse\n" if(exists $anno{$queryTaxID});
	}
}

# Fix some missing verification results
open IN,"$outdir/validated.fna" or die $!;
while(<IN>){
	chomp;
	next unless(/^>/);
	my $id = $_;
	$id =~ s/^>//;
	my $queryTaxID = (split /\|/,$id)[-1];
	print VO "$id\t$anno{$queryTaxID}\tfalse\n" if(!exists $score{$id} && exists $anno{$queryTaxID});
}
close IN;
close VO;

sub help{
print "
	Usage: perl $0

	--report             kraken's report file
	--inMark             enter a list of sequences id to be validated.
	--db                 path of the database to annotate
	--fq                 clean data's fastq file, separated with comma, example: --fq fq1,fq2
	--krakenout          intermediate file of kraken analysis, the file's suffix is out
	--outdir             output path. [./]
	--version|v          print version information.

";
exit(0);
}
