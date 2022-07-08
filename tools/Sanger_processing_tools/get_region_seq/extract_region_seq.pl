use strict;
#use warn;
my $fa_dir="/results/zhongkejiying/databases/ucsc/hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips";
my $refGene="/results/zhongkejiying/databases/ucsc/hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt";
my ($chr,$pos1,$pos2)=split(/:/,$ARGV[1]);
my $line = "";
my $reaction_maxLen = 1500;
my $least_dist2Start = 500;
my %h_info = ();
my $last_project = "";
while( $line = <STDIN> ){
	chomp $line;
	next if( $line =~ /^#/ );
	if( $line =~ /^>(.*)$/ ){
		$last_project = $1;
		next;
	}
	my @arr = split(/\t/,$line);
	my ( $chr,$pos1,$pos2 ) = split(/:/,$arr[0]);

	if( $chr =~ /M/i ){
		print STDERR "ERR:chrM is not supported.\n$line\n";
	}
	if( $chr =~ /[xy]/ ){
		print STDERR "ERR:chrx,chry should be input sa chrX,chrY.\n$line\n";
	}

	$h_info{$last_project}{$chr}{$pos1}{$pos2} = $arr[1];
}


foreach my $project( sort {$a cmp $b} keys(%h_info) ){

	print ">$project\n";
	my %h_regions = ();
	foreach my $chr( sort {$a cmp $b} keys(%{$h_info{$project}}) ){

		#waste setting to get the data all output
		$h_info{$project}{$chr}{10000000000000}{100000000000000}++;

		#get reference sequence
		my  $chr_fa = "";
		open REF,"$fa_dir/$chr.fa" or die $!;
		$/ = ">";
		while( $line = <REF> ){
			if( $line =~ s/^(\S+)\s// ){
				my $t_chr = $1;
				if( $chr eq $t_chr ){
					$line =~ s/\s//g;
					$chr_fa = $line;
					last;
				}
			}
		}
		$/ = "\n";
		close REF;

		my $last_region = -1000000;
		foreach my $pos1( sort {$a <=> $b} keys(%{$h_info{$project}{$chr}}) ){
			foreach my $pos2( sort {$a <=> $b} keys(%{$h_info{$project}{$chr}{$pos1}}) ){
				if( ( $last_region + $reaction_maxLen - $least_dist2Start ) < $pos1 ){
					if( $last_region > 0 ){
						my $region_len = $h_regions{ $last_region }{"last_pos"}-$last_region+$least_dist2Start+1;
						my $mut_count = 0;
						foreach my $t_pos1( sort {$a <=> $b} keys( %{$h_regions{ $last_region }{"mut"}} ) ){
							foreach my $t_pos2( sort {$a <=> $b} keys( %{$h_regions{ $last_region }{"mut"}{$t_pos1}} ) ){
								$mut_count++;
							}
						}
						print "Region_info:\tmut_count:$mut_count\t$chr:$region_len:$last_region:".($h_regions{ $last_region }{"last_pos"}+$least_dist2Start)."\n";

						#mark mutations
						my $count = 0;
						my $t_chr_fa = $chr_fa;
						foreach my $t_pos1( sort {$a <=> $b} keys( %{$h_regions{ $last_region }{"mut"}} ) ){
							foreach my $t_pos2( sort {$a <=> $b} keys( %{$h_regions{ $last_region }{"mut"}{$t_pos1}} ) ){
								$count++;
								my $t_len = $t_pos2 - $t_pos1+1;
								$t_chr_fa = substr($t_chr_fa,0,$t_pos1-1)."X"x$t_len.substr($chr_fa,$t_pos2);
								my $ref = substr($chr_fa,$t_pos1,$t_len);
								print "mutation$count:\tref:$chr:$t_pos1:$t_pos2:$ref\t".$h_regions{ $last_region }{"mut"}{$t_pos1}{$t_pos2}."\n";
							}
						}
						my $out_str = substr($t_chr_fa,$last_region,$region_len);

						#print sequence by 50 characters per line
						my $count_per_line = 100;
						while($out_str =~ s/(.{$count_per_line})//){
							my $t_out = $1;
							print $t_out."\n";
						}
						print $out_str."\n\n";

					}
					$last_region = $pos1 - $least_dist2Start;
					$h_regions{ $last_region }{"mut"}{$pos1}{$pos2} = $h_info{$project}{$chr}{$pos1}{$pos2};
					$h_regions{ $last_region }{"last_pos"} = $pos2;
				}
				else{
					if( $h_regions{ $last_region }{"last_pos"} < $pos2 ){
						$h_regions{ $last_region }{"last_pos"} = $pos2;
					}
					$h_regions{ $last_region }{"mut"}{$pos1}{$pos2}= $h_info{$project}{$chr}{$pos1}{$pos2};
				}
			}
		}
	}
}
