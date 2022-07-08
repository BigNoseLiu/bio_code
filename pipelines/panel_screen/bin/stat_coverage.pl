my %h_stat_total = ();
my %h_stat = ();
my %h_types= ();
my %h_out_panel = ();
my %h_target_panel = ();
while( $file = <STDIN> ){
	chomp $file;
	if( $file =~ /\/([^\/]+)\.([^\.]+)_cover_(.*)$/ ){
		my $panel = $1;
		my $type = $2;
		my $out_panel = $3;
		$h_types{$type}++;
		$h_target_panel{$panel}++;
		$h_out_panel{$out_panel}++;
		open IN, $file or die $!;
		while( $line = <IN> ){
			chomp $line;
			my @arr = split(/\t/,$line);
			$h_stat{$arr[3]}{'out'}{$out_panel}{$type}{'total'} += $arr[7];
			$h_stat{$arr[3]}{'out'}{$out_panel}{$type}{'cover'} += $arr[6];
			$h_stat{$arr[3]}{'target'}{$panel}++;
			$h_stat_total{$panel}{$out_panel}{$type}{'total'} += $arr[7];
			$h_stat_total{$panel}{$out_panel}{$type}{'cover'} += $arr[6];
		}
		close IN;
	}
}

print "#gene";
foreach $out_panel( sort {$a cmp $b} keys(%h_out_panel) ){
	foreach $type( sort {$a cmp $b} keys(%h_types) ){
		print "\t$out_panel\_$type";
	}
}
foreach $panel( sort {$a cmp $b} keys(%h_target_panel) ){
	print "\t$panel";
}
print "\n";



foreach my $gene( sort {$a cmp $b}  keys(%h_stat)){
	print $gene;
	foreach $out_panel( sort {$a cmp $b} keys(%h_out_panel) ){
		foreach $type( sort {$a cmp $b} keys(%h_types) ){
			if( defined($h_stat{$gene}{'out'}{$out_panel}{$type}{'total'}) ){
				my $cover = $h_stat{$gene}{'out'}{$out_panel}{$type}{'cover'}/$h_stat{$gene}{'out'}{$out_panel}{$type}{'total'};
				print "\t$cover";
				if( $cover < 0.8 ){
					foreach $panel( sort {$a cmp $b} keys(%h_target_panel) ){
						$h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}{$gene} = $cover;
					}
				}
			}
			else{
				print "\t-";
			}
		}
	}
	foreach $panel( sort {$a cmp $b} keys(%h_target_panel) ){
		if( defined($h_stat{$gene}{'target'}{$panel}) ){
			print "\t1";
		}
		else{
			print "\t-";
		}
	}
	print "\n";
}


print STDERR "#target_panel";
foreach $out_panel( sort {$a cmp $b} keys(%h_out_panel) ){
	foreach $type( sort {$a cmp $b} keys(%h_types) ){
		print STDERR "\t$out_panel\_$type";
	}
}
print STDERR "\n";


foreach $panel( sort {$a cmp $b} keys(%h_target_panel) ){
	print STDERR $panel;
	foreach $out_panel( sort {$a cmp $b} keys(%h_out_panel) ){
		foreach $type( sort {$a cmp $b} keys(%h_types) ){
			if( defined($h_stat_total{$panel}{$out_panel}{$type}{'total'}) ){
				my $cover = $h_stat_total{$panel}{$out_panel}{$type}{'cover'}/$h_stat_total{$panel}{$out_panel}{$type}{'total'};
				print STDERR "\t$cover";
				if( defined($h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}) ){
					my @arr_lack_genes =  sort {$h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}{$b} <=> $h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}{$a}} keys(%{$h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}});
					print STDERR "(".scalar(@arr_lack_genes).")";
					foreach my $gene( @arr_lack_genes ){# sort {$h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}{$b} <=> $h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}{$a}} keys(%{$h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}}) ){
						print STDERR "|$gene:$h_stat_total{$panel}{$out_panel}{$type}{'lack_genes'}{$gene}";
					}
				}
			}
			else{
				print STDERR "\t-";
			}
		}
	}
	print STDERR "\n";
}
