use strict;
my $line = "";
my $indent_char = " ";
my $current_indentChar = "";
while( $line = <STDIN> ){
	chomp $line;
	$line =~ s/></>\n</g;
	my @arr = split(/\n/,$line);
	foreach my $t_line ( @arr ){
		if( $t_line =~ /^<\// ){
			$current_indentChar =~ s/$indent_char$//;
		}
		print $current_indentChar.$t_line."\n";
		$current_indentChar .= $indent_char if( $t_line !~ /<\// && $t_line !~ /\/>$/ );
	}
}
