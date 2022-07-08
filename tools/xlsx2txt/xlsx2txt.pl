#!/usr/bin/env perl
use strict;
use warnings;
use Spreadsheet::Read;
use Encode;
use File::Basename;
use File::Spec::Functions qw(rel2abs);

@ARGV==2 || die "Usage: perl $0 <input.xlsx> <sheet_dir>\n";
#my $type = ($ARGV[0] =~/xlsx$/)? "xlsx" : "xls";
my $out_dir = $ARGV[1];
`mkdir -p $out_dir`;
my $outname=basename ($ARGV[0]);
&parse_xlsx($ARGV[0]);

sub parse_xlsx {
	my @file_lines=();
	my $file_path = shift;
	my $file_type = shift;
	my $workbook = ReadData($file_path,cells => 0, dtfmt => "yyyy-mm-dd", parser => "xlsx");
	if(defined $workbook->[0]{'error'}){
		print "Error occurred while processing $file_path:".
			$workbook->[0]{'error'}."\n";
		exit(-1);
	}
	my $sheet_count = $workbook->[0]{"sheets"} or die "No sheets in $file_path/n";
	for my $sheet_index (1 .. $sheet_count){
		my $worksheet = $workbook->[$sheet_index];
		my $sheet_name = $worksheet->{'label'};
		open FO,">$out_dir/$outname.$sheet_name.txt";
		my $max_rows = $worksheet->{'maxrow'};
		my $max_cols = $worksheet->{'maxcol'};
		for my $row_num (1..($max_rows)){
			my @tmp = ();
			for my $col_num (1..($max_cols)){
				my $cell = $worksheet->{'cell'}[$col_num][$row_num];
				if (defined $cell){
					$cell =~ s/\n|\r/^/g;
					$cell =~ s/^\s+//g;
					$cell =~ s/\s+$//g;
					if ($cell eq ''){
						push @tmp,'';
					}else{
						push @tmp,$cell;
					}
				}else{
					push @tmp,'';
				}
			}
			print FO (join "\t",@tmp)."\n";
		}
		close FO;
	}
}
