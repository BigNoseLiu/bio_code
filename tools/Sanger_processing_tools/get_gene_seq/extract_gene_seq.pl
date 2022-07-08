use strict;
#use warn;
my $fa_dir="/results/zhongkejiying/databases/ucsc/hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips";
my $refGene="/results/zhongkejiying/databases/ucsc/hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt";
my ($chr,$pos1,$pos2)=split(/:/,$ARGV[1]);
my $line = "";
my %h_info = ();
my %h_nm_id = ();
my $last_project = "";
while( $line = <STDIN> ){
	chomp $line;
	next if( $line =~ /^#/ );
	if( $line =~ /^>(.*)$/ ){
		$last_project = $1;
		next;
	}
	my @arr = split(/\t/,$line);
	my ( $chr,$pos1,$pos2,$nm_id ) = split(/:/,$arr[0]);

	if( $chr =~ /M/i ){
		print STDERR "ERR:chrM is not supported.\n$line\n";
	}
	if( $chr =~ /[xy]/ ){
		print STDERR "ERR:chrx,chry should be input sa chrX,chrY.\n$line\n";
	}

	$h_info{$last_project}{$chr}{$nm_id}{$pos1}{$pos2} = $arr[1];
	$h_nm_id{$nm_id} = 0;
}

#get ref_gene info
open GENE,$refGene or die $!;
my @lines = <GENE>;
foreach $line ( @lines ){
	chomp $line;
	my @arr = split(/\t/,$line);
	if( defined( $h_nm_id{$arr[1]} ) ){
		$h_nm_id{$arr[1]} = $line;
	}
}
close GENE;

foreach my $project( sort {$a cmp $b} keys(%h_info) ){
	print ">Project_name:\t$project\n";

	foreach my $chr( sort {$a cmp $b} keys(%{$h_info{$project}}) ){
		my  $chr_fa = "";
		#get reference sequence
		open REF,"$fa_dir/$chr.fa" or die $!;
		$/ = ">";
		while( $line = <REF> ){
			if( $line =~ s/^(\S+)\s// ){
				my $t_chr = $1;
				if( $chr eq $t_chr ){
					$line =~ s/\s//g;
					$chr_fa = $line;
					last;
					my $out_str = substr($line,$pos1-1,$pos2-$pos1+1);
				}
			}
		}
		$/ = "\n";
		close REF;

		#each nm_id is a print unit
		foreach my $nm_id( sort {$a cmp $b} keys(%{$h_info{$project}{$chr}}) ){
			if( $h_nm_id{$nm_id} == 0 ){
				print STDERR "$nm_id not found\n";
				next;
			}

			#get nm_id position info
			my @t_arr = split(/\t/,$h_nm_id{$nm_id});
			my @pos1_s = split(/,/,$t_arr[9]);
			my @pos2_s = split(/,/,$t_arr[10]);
			print "gene_info:\t".join(":",@t_arr[12,1,3,2,4,5])."\n";

			#mark mutation position
			my $mut_count = 0;
			foreach my $pos1( sort {$a <=> $b} keys(%{$h_info{$project}{$chr}{$nm_id}}) ){
				foreach my $pos2( sort {$a <=> $b} keys(%{$h_info{$project}{$chr}{$nm_id}{$pos1}}) ){
					$mut_count++;
					my $len = $pos2-$pos1+1;
					my $info = $h_info{$project}{$chr}{$nm_id}{$pos1}{$pos2};
					print "mutation$mut_count:\tref:$chr:$pos1:$pos2:".substr($chr_fa,$pos1-1,$pos2-$pos1+1)."\t".$h_info{$project}{$chr}{$nm_id}{$pos1}{$pos2}."\n";
					$chr_fa = substr($chr_fa,0,$pos1-1)."X"x$len.substr($chr_fa,$pos2);
				}
			}

			#get the sequence of exons & introns
			my $str_out = "";
			foreach( my $i=0;$i<scalar @pos1_s;$i++ ){
				if( $pos1_s[$i] !~/\d/ ){
					last;
				}
				if( $i>0 ){
					#set introns sequence as lower charater
					$str_out .= lc(substr($chr_fa,$pos2_s[$i-1],$pos1_s[$i]-$pos2_s[$i-1]));
				}
				#set exons sequence as upper charater
				$str_out .= uc(substr($chr_fa,$pos1_s[$i],$pos2_s[$i]-$pos1_s[$i]));
			}

			#print sequence by 50 characters per line
			my $count_per_line = 50;
			while($str_out =~ s/(.{$count_per_line})//){
				my $t_out = $1;
				print $t_out."\n";
			}
			print $str_out."\n\n";
		}
	}
	print "\n\n\n";
}
