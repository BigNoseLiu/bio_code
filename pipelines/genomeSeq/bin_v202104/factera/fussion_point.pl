#input perl $0 sample.markDup.sorted.bam 14:65239442 16:8424213 out_file
my ( $bam_file,$break1,$break2,$out_file ) = @ARGV;
my ($break1_chr,$break1_pos) = split(/:/,$break1);
my ($break2_chr,$break2_pos) = split(/:/,$break2);
my $break1_region = "$break1_chr:".($break1_pos-300)."-".($break1_pos+300);
my $break2_region = "$break2_chr:".($break2_pos-300)."-".($break2_pos+300);
`samtools view -F 4 $bam_file $break1_region >$out_file\.temp`;
`samtools view -F 4 $bam_file $break2_region >>$out_file\.temp`;
`cat $out_file\.temp|sort|uniq >$out_file\.sam`;
`rm -rf $out_file\.temp`;

my %h_stat = ();
my $line;
open IN,"$out_file\.sam" or die $!;
while($line = <IN>){
	chomp $line;
	my @arr = split( /\t/, $line );
	my $sam_flag = sprintf("%016b",$arr[1]);
	my $strand = "+";
	$strand = "-" if($sam_flag =~ /1....$/);
	$h_stat{$arr[0]}{$arr[2]."\t".$arr[3]."\t".$strand}++;
}
close IN;

open OUT,">$out_file" or die $!;
foreach my $read( sort {$a cmp $b} keys(%h_stat) ){
	my @positions = sort {$a cmp $b} keys(%{$h_stat{$read}});
	my $pos_count = scalar@positions;
	if( $pos_count < 2 ){
		next;
	}
	for( my$i=0;$i<(scalar@positions-1);$i++ ){
		my $p1 = $positions[$i];
		for( my$j=$i+1;$j<scalar@positions;$j++ ){
			my $p2 = $positions[$j];
			print OUT "$read\t$pos_count\t$p1\t$p2\n";
		}
	}
}
close OUT;
