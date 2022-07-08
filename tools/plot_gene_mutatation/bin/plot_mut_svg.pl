#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(rmtree);
use Cwd qw(abs_path);
use Data::Dumper;

sub usageWithColor {
	print color "green";#change the text color
	print <<USAGE;
Description:
        Used to plot the mutations.
Usage:  
        perl $0 [options]
        Options:
	     -help : reveal help info
	     -out  <str>  : output file path
        Example:
             perl $0 -out  output_file
Author & Contact:
	Mingming Liu
        mingming_liu\@shengtingroup.com
Last updated:
        2013-10-16
USAGE
	print color "reset";#change back the text color
}
sub err_print{
	print STDERR color "red";#change the text color
	print STDERR $_[0];
	print STDERR color "reset";#change back the text color
}

my $GRCh_version = "37";
my ($help,$infile,$outfile,$only_chr);
my $line;
my ( $total_width,$total_height,$total_height2 ) = (800,400,1000);
my $chr_len_max = 250000000;
my $arial_10_scale = (330/50);
#my $chinese_font = "Times New Roman";
#my $chinese_font = "Arial";
#my $chinese_font = "Arial";
#my $chinese_font = "华文楷体";
my $chinese_font = "楷体";

my %h_cytoband_color1 = (
	"acen"=>"rgb(255,106,106)",
	"gneg"=>"rgb(238,213,183)",
	"gpos25"=>"rgb(192,255,62)",
	"gpos50"=>"rgb(179,238,58)",
	"gpos75"=>"rgb(154,205,50)",
	"gpos100"=>"rgb(105,139,34)",
	"gvar"=>"rgb(238,230,133)",
	"stalk"=>"rgb(192,255,62)",
);



my %h_cytoband_color2 = (
	"acen"=>"rgb(139,101,8)",
	"gneg"=>"rgb(152,245,255)",
	"gpos25"=>"rgb(132,112,255)",
	"gpos50"=>"rgb(123,104,238)",
	"gpos75"=>"rgb(106,90,205)",
	"gpos100"=>"rgb(72,61,139)",
	"gvar"=>"rgb(0,0,0)",
	"stalk"=>"rgb(192,255,62)",
);
my %h_cytoband_color = %h_cytoband_color1;
GetOptions(
	"help"=>\$help,
	"in=s"=>\$infile,
	"out=s"=>\$outfile,
	"GRCh=s"=>\$GRCh_version,
	"only"=>\$only_chr,
);
if (defined $help ) {
	&usageWithColor();
	exit 0;
}
#starting
print "Start $0 time:".localtime(time)."\n";
#read databases
my ($cytoband_file,$gff_file) = ("","");
open CONF,"../config/database.list" or die $!;
while( $line = <CONF> ){
	chomp $line;
	my @arr = split(/=/,$line);
	if( $arr[0] =~ /^\s*GRCh$GRCh_version\_gff\s*$/ ){
		$gff_file = $arr[1];
	}
	elsif( $arr[0] =~ /^\s*GRCh$GRCh_version\_cytoBand\s*$/ ){
		$cytoband_file = $arr[1];
	}
}
close CONF;
#read cytoband
my %h_cyto;
my %h_chr_attr = ();
open CYTO,"gunzip -c $cytoband_file|" or die $!;
while( $line = <CYTO> ){
	chomp $line;
	my @arr = split(/\s+/,$line);
	$h_cyto{ $arr[0] }{$arr[1]}{"end"} = $arr[2];
	$h_cyto{ $arr[0] }{$arr[1]}{"pq"} = $arr[3];
	$h_cyto{ $arr[0] }{$arr[1]}{"color"} = $arr[4];
	if( $arr[3] =~ /^p/ ){
		$h_chr_attr{$arr[0]}{"p_start"} = $arr[1] if( !defined($h_chr_attr{$arr[0]}{"p_start"}) );
		$h_chr_attr{$arr[0]}{"p_end"} = $arr[1];
	}
	elsif( $arr[3] =~ /^q/ ){
		$h_chr_attr{$arr[0]}{"q_start"} = $arr[1] if( !defined($h_chr_attr{$arr[0]}{"q_start"}) );
		$h_chr_attr{$arr[0]}{"q_end"} = $arr[1];
	}
	$h_chr_attr{$arr[0]}{"len"} = $arr[2];
}
close CYTO;

#read gff
open GFF,"gunzip -c $gff_file|" or die $!;
my %h_id2gene = ();
my %h_id2mRNA = ();
my %h_mRNA_attr = ();
my %h_exons = ();
while( $line = <GFF> ){
	chomp $line;
	next if( $line =~ /^#/ );
	my @arr = split(/\t/,$line);
	if( $line =~ /\tgene\t.*\tID=([^;]+);Name=([^;]+);/ ){
		$h_id2gene{ $1 }{"name"} = $2;
	}
	elsif( $arr[2] eq "mRNA" || $arr[2] eq "transcript"  || $arr[2] eq "ncRNA"){
		$line =~ /ID=([^;]+);Name=([^\.]+)\.\d+;Parent=([^;]+);/ ;
		my ( $id,$name,$parent ) = ($1,$2,$3);
		$h_id2mRNA{ $id }{"name"} =  $name;
		$h_mRNA_attr{$name}{"gene"} = $h_id2gene{ $parent }{"name"};
		$h_mRNA_attr{$name}{"start"} = $arr[3];
		$h_mRNA_attr{$name}{"end"} = $arr[4];
		$h_mRNA_attr{$name}{"len"} = ($arr[4] - $arr[3]);
		$h_mRNA_attr{$name}{"strand"} = $arr[6];
	}
	elsif( $line =~ /\texon\t.*\tID=([^;]+);Parent=([^;]+);/ ){
		my ( $id,$parent ) = ($1,$2);
		if( defined( $h_id2mRNA{ $parent }{"name"} ) ){
			my $rna_name = $h_id2mRNA{ $parent }{"name"};
			$h_exons{$rna_name}{$arr[3]} = $arr[4];
		}
	}
}
close GFF;
&plot("deafness.1.svg","chr13","GJB2","NM_004004","rs80338942","NM_004004.5:c.167delT,NC_000013.10:g.20763554del","Deletion","Leu56Argfs",20763554);
#ending
print "End   $0 time:".localtime(time)."\n";
sub plot{
	my ( $svg_file,$chr,$gene_name,$nm_id,$rs_id,$hgvs,$mut_type,$protein_change,$mut_pos ) = @_;
	#save STDOUT handle
	open STD,">&STDOUT" or die $!;

	#print svg file
	open STDOUT,">$svg_file" or die $!;
	#&svg_begin_with_width_height($total_width,$total_height);
	&svg_begin();
	#&svg_rect( 5,5,20,20,$total_width-5,$total_height-5,"white","red",1);
	&plot_chr($chr);
	&plot_gene($nm_id,$chr,$gene_name);
	&plot_mutation($rs_id,"NA",$nm_id,$mut_pos,$mut_type,$protein_change,$hgvs);
	&title("突变: $rs_id");
	&svg_end();
	close STDOUT;

	#get STDOUT back
	open STDOUT,">&STD" or die $!;
	close STD;

}
sub svg_begin_with_width_height{
	my ($width,$height) = @_;
	print <<SVG;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="$width" height="$height" version="1.1" xmlns="http://www.w3.org/2000/svg">
SVG
}
sub svg_begin{
	#my ($width,$height) = @_;
	print <<SVG;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="100%" height="100%" version="1.1" xmlns="http://www.w3.org/2000/svg">
SVG
}
#<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "liumingming">
sub svg_end{
	print "</svg>";
}

sub svg_rect_with_opacity{
	my ($x,$y,$rx,$ry,$width,$height,$fill,$stroke_col,$stroke_width,$opacity) = @_;
	print <<SVG;
<rect x="$x" y="$y" rx="$rx" ry="$ry" width="$width" height="$height" style="fill:$fill;stroke:$stroke_col;stroke-width:$stroke_width;opacity:$opacity"/>
SVG
}

sub svg_rect{
	my ($x,$y,$rx,$ry,$width,$height,$fill,$stroke_col,$stroke_width) = @_;
	print <<SVG;
<rect x="$x" y="$y" rx="$rx" ry="$ry" width="$width" height="$height" style="fill:$fill;stroke:$stroke_col;stroke-width:$stroke_width;opacity:1"/>
SVG
}
sub svg_circle{
	my ($cx,$cy,$r,$stroke_col,$stroke_width,$fill) = @_;
	print <<SVG;
<circle cx="$cx" cy="$cy" r="$r" stroke="$stroke_col" stroke-width="$stroke_width" fill="$fill"/>
SVG
}
sub svg_ellipse{
	my ($cx,$cy,$rx,$ry,$fill,$stroke_col,$stroke_width,$opacity) = @_;
	print <<SVG;
<ellipse cx="$cx" cy="$cy" rx="$rx" ry="$ry" style="fill:$fill;stroke:$stroke_col;stroke-width:$stroke_width;opacity:$opacity"/>
SVG
}


sub svg_polygon{
	my ($points,$fill,$stroke_col,$stroke_width,$opacity) = @_;
	print <<SVG;
<polygon points="$points" style="fill:$fill;stroke:$stroke_col;stroke-width:$stroke_width;opacity:$opacity"/>
SVG
}



sub svg_polyline{
	my ($points,$stroke_col,$stroke_width) = @_;
	print <<SVG;
<polyline points="$points" style="fill:white;stroke:$stroke_col;stroke-width:$stroke_width"/>
SVG
}

sub svg_text_rotate{
	my ( $x,$y,$color,$text,$font,$font_size,$rotate,$weight) = @_;
	print <<SVG;
<text x="$x" y="$y" fill="$color" font-family="$font" font-size="$font_size" transform="rotate($rotate)" style="font-weight:$weight;">$text</text>
SVG
#<text x="$x" y="$y" fill="$color" font-family="$font" font-size="$font_size" transform="rotate($rotate)">$text</text>
}


sub svg_text{
	my ( $x,$y,$color,$text,$font_size ) = @_;
	print <<SVG;
<text x="$x" y="$y" fill="$color" font-family="Arial" font-size="$font_size">$text</text>
SVG
}

sub svg_line{
	my ($x1,$y1,$x2,$y2,$stroke_col,$stroke_width) = @_;
	print <<SVG;
<line x1="$x1" y1="$y1" x2="$x2" y2="$y2" style="stroke:$stroke_col;stroke-width:$stroke_width"/>
SVG
}
sub title{
	my ( $title_text ) = @_;
	my $scale = 3.5;
	&svg_text_rotate($total_width/15,35,"dark",$title_text,"楷体",30,"0 0,0","常规");
	#&svg_text_rotate($total_width/15,35,"dark",$title_text,$chinese_font,30,"0 0,0","bold");
	#&svg_text_rotate($total_width/2-length($title_text)*$scale,35,"dark",$title_text,$chinese_font,24,"0 0,0","bold");
}
sub per_thousand_sign{
	my ($num_in) = @_;
	my $number = "$num_in";
	my @out;
	while( $number =~ s/(...)$// ){
		push @out,$1;
	}
	push @out,$number if( $number =~ /./ );
	my @out_num;
	for(my $i=scalar @out;$i>0;$i--){
		push @out_num,$out[$i-1];
	}
	return join(",",@out_num);
}
sub get_gene_pos{
	my ($nm_id,$pos1,$pos2) = @_;
	my $x_start = (3/10)*$total_width;
	my $y_start = (1/30+2/5)*$total_height;
	my $scale = (9*$total_width/10-$x_start)/$h_mRNA_attr{$nm_id}{"len"};
	my $x = $x_start + ($pos1-$h_mRNA_attr{$nm_id}{"start"})*$scale;
	my $height = $total_height/20;
	my $width= $scale*($pos2-$pos1);
	return ( $x,$y_start,$width,$height );
}


sub get_chr_pos{
	my ($chr,$pos1,$pos2) = @_;
	my $x_start = $total_width/10;
	my $y_start = $total_height/5;
	my $scale = (19*$total_height2/20-$y_start)/$h_chr_attr{$chr}{"len"};
	#my $scale = (19*$total_height/20-$y_start)/$h_chr_attr{$chr}{"len"};
	my $y = $y_start + $pos1*$scale;
	my $width = $total_width/20;
	my $height = $scale*($pos2-$pos1);
			#&svg_rect( $x,$last_y,$rx,$ry,$width,$height,$fill_col,"white",0);
	return ( $x_start,$y,$width,$height );
}

sub plot_mutation{
	my ($rs_id,$mut_name,$nm_id,$pos,$mut_type,$protein_change,$hgvss) = @_;
	my ($gene_start_x,$gene_start_y,$gene_width,$gene_height) = &get_gene_pos($nm_id,$h_mRNA_attr{$nm_id}{"start"},$h_mRNA_attr{$nm_id}{"end"});
	my ($mut_start_x,$mut_start_y,$mut_width,$mut_height) = ( $gene_start_x+50,$gene_start_y+$gene_height+40,350,130 );
	&svg_rect_with_opacity( $mut_start_x,$mut_start_y,15,15,$mut_width,$mut_height,"red","black",0,0.4);
	my ($title_text_x ,$title_text_y) = ($mut_start_x+10,$mut_start_y+25);
	&svg_text_rotate($title_text_x,$title_text_y,"dark","突变: ",$chinese_font,18,"0 0,0","bold");
	my $mut_font_size = 14;
	&svg_text_rotate($title_text_x+50,$title_text_y,"dark","$rs_id","Arial",$mut_font_size,"0 0,0","bold") if($rs_id ne "-");
	&svg_text_rotate($title_text_x+30,$title_text_y+25,"dark","突变类型：","宋体",$mut_font_size,"0 0,0","bold");
	&svg_text_rotate($title_text_x+110,$title_text_y+25,"dark",$mut_type,"宋体",$mut_font_size,"0 0,0","normal");
	&svg_text_rotate($title_text_x+30,$title_text_y+45,"dark","蛋白影响：","宋体",$mut_font_size,"0 0,0","bold");
	&svg_text_rotate($title_text_x+110,$title_text_y+45,"dark",$protein_change,"宋体",$mut_font_size,"0 0,0","normal");
	&svg_text_rotate($title_text_x+30,$title_text_y+65,"dark","HGVS命名：","宋体",$mut_font_size,"0 0,0","bold");
	my $hgvs_y = $title_text_y+65;
	foreach my $hgvs( split(/,/,$hgvss) ){
		&svg_text_rotate($title_text_x+110,$hgvs_y,"dark",$hgvs,"宋体",$mut_font_size,"0 0,0","normal");
		$hgvs_y += 15;
	}
	my ($x,$y,$width,$height) = &get_gene_pos($nm_id,$pos,$pos);
	my ($x1,$y1) = ($x-8,$y+$height+40);
	$x1 = ($mut_start_x+10) if( $x1 < ($mut_start_x+20) );
	$x1 = ($mut_start_x+$mut_width-10) if( $x1 > ($mut_start_x+$mut_width-20) );
	my ($x2,$y2) = ($x,$y+$height);
	my ($x3,$y3) = ($x1+10,$y1);
	&svg_polygon("$x1,$y1 $x2,$y2 $x3,$y3","red","blue",0,0.4);
	&svg_ellipse($x2,$y+$height/2,2,$height-2,"red","red",0,0.5);
	#&svg_rect_with_opacity( $x-5,$y,3,3,20,$height,"red","dark",3,0.5);
	my $text_x = $x1;
	my $text = "Mutation: $mut_name";
	my $text_x_limit = $gene_start_x+$gene_width-$arial_10_scale*length($text);
	$text_x =  $text_x_limit if( $text_x >$text_x_limit );
	#&svg_rect_with_opacity( $text_x-5,$y1,3,3,$arial_10_scale*length($text),20,"red","black",0,0.4);
	#&svg_text($text_x,$y1+15,"black",$text,10);
}
sub plot_gene{
	my ($nm_id,$chr,$gene_name) = @_;
	my $stroke_col = "black";
	my $stroke_width = 2;
	my ($rx,$ry) = (0,0);
	#mark gene position on chr
	my ($chr_x,$chr_y,$chr_width,$chr_height) = &get_chr_pos($chr,$h_mRNA_attr{$nm_id}{"start"},$h_mRNA_attr{$nm_id}{"end"});
	&svg_rect_with_opacity( $chr_x-10,$chr_y-5,3,3,$chr_width+20,$chr_height+10,"blue",$stroke_col,0,0.5);
	#&svg_rect_with_opacity( $chr_x-10,$chr_y-5,3,3,$chr_width+20,$chr_height+10,"rgb(255,165,0)",$stroke_col,0,0.5);
	
	#plot gene 
	my ($x,$y,$width,$height) = &get_gene_pos($nm_id,$h_mRNA_attr{$nm_id}{"start"},$h_mRNA_attr{$nm_id}{"end"});
	&svg_rect_with_opacity( $x-30,$y-$height-90,10,10,$width+70,140,"blue",$stroke_col,0,0.2);
	my ($x1,$y1) = ($x-30,$y-$height-70);
	my ($x2,$y2) = ($chr_x+$chr_width+5,$chr_y);
	my ($x3,$y3) = ($chr_x+$chr_width+5,$chr_y+$chr_height);
	my ($x4,$y4) = ($x-30,$y-$height-40);
	&svg_polygon("$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4","blue","blue",0,0.2);
	#plot gene pos
	&svg_rect( $x,$y-$height,$rx,$ry,$width,1,"black",$stroke_col,0);
	#gene start pos
	&svg_rect( $x,$y-3*$height/2,$rx,$ry,1,$height/2,"black",$stroke_col,0);
	&svg_text($x-5,$y-3*$height/2-8,"black","From: ".&per_thousand_sign($h_mRNA_attr{$nm_id}{"start"}),14);
	#gene end pos
	&svg_rect( $x+$width-1,$y-3*$height/2,$rx,$ry,1,$height/2,"black",$stroke_col,0);
	my $to_pos = &per_thousand_sign($h_mRNA_attr{$nm_id}{"end"});
	my $arial_8_scale = (300/50);
	&svg_text($x+$width-$arial_8_scale*length($to_pos)-10,$y-3*$height/2-8,"black","To: $to_pos",14);
	#title
	#&svg_text($x+$width/2-20,$y-3*$height/2-28,"rgb(139,69,19)",$nm_id,20);
	my ($text_gene,$text_x,$text_y) = ("基因: ".$gene_name,$x+$width/2-60,$y-3*$height/2-50);
	&svg_text_rotate($text_x,$text_y,"dark",$text_gene,$chinese_font,18,"0 0,0","bold");

	#plot mark area
	#&svg_line($chr_x+$chr_width+10,$chr_y+$chr_height/2,$x,$y-3*$height/2-60,"green",2);
	#&svg_line($chr_x+$chr_width+10,$chr_y+$chr_height/2,$x-15,$y+$height,"green",2);

	#plot exons 
	my @positions = sort {$a <=> $b} keys(%{$h_exons{$nm_id}});
	my ($last_x) = (-1);
	foreach my $pos1(@positions){
		($x,$y,$width,$height) = &get_gene_pos($nm_id,$pos1,$h_exons{$nm_id}{$pos1});
		&svg_rect( $x,$y,$rx,$ry,$width,$height,"green",$stroke_col,0);
		if( $last_x > 0 ){
			my $y2 = $y+$height/2;
			my $x2 = ($x+$last_x)/2;
			my $y3 = $y2-2;
			&svg_polyline("$x,$y2 $x2,$y3 $last_x,$y2","red",1);
			#&svg_polyline("$x,$y2 $x2,$y $last_x,$y2","red",1);
		}
		$last_x = $x+$width;
	}
}


sub plot_chr{
	my ($chr) = @_;
	#title
	my ($x,$y,$width,$height) = &get_chr_pos($chr,0,$h_chr_attr{$chr}{"len"});
	my $text_chr = $chr;
	$text_chr =~ s/^chr(.*)$/$1 染色体/;
	$text_chr =~ s/^(\d+)([^\d])/$1 号/;
	#$text_chr =~ s/chr/Chromosome /;
	my ($text_x,$text_y) = ($x-25,300);
	#my ($text_x,$text_y) = ($x-20,$y+$height/2+25);
	&svg_text_rotate($text_x,$text_y,"dark",$text_chr,$chinese_font,22,"270 $text_x,$text_y","bold");
	#&svg_text(15,$y-15,"green",$text_chr,22);

	my @arr_last = ();
	my $stroke_col = "black";
	my $stroke_width = 2;
	my ($rx,$ry) = (0,0);
	my ($line_start_x,$line_start_y);
	my $ellipse_max_height = 10;
	my @chips_pos = ();
	foreach my $attr( ("p_start","p_end","q_start","q_end") ){
		my $pos1 = $h_chr_attr{$chr}{ $attr };
		($x,$y,$width,$height) = &get_chr_pos($chr,$pos1,$h_cyto{$chr}{$pos1}{"end"});
		my $fill_col = $h_cytoband_color{ $h_cyto{$chr}{$pos1}{"color"} };
		#if( $h_chr_attr{$chr}{"p_start"} == $pos1 || $h_chr_attr{$chr}{"q_start"} == $pos1 ){
		my $rect_height = 0;
		if( $height > $ellipse_max_height ){
			$rect_height = $height - $ellipse_max_height;
			$height = $ellipse_max_height;
		}
		if( $attr eq "p_start" || $attr eq "q_start" ){
			&svg_ellipse($x+$width/2,$y+$height,$width/2,$height,$fill_col,$stroke_col,1,1);
			&svg_line($x,$y+$height,$x,$y+$height+$rect_height,$stroke_col,$stroke_width);
			&svg_line($x+$width,$y+$height,$x+$width,$y+$height+$rect_height,$stroke_col,$stroke_width);
			&svg_rect( $x,$y+$height,$rx,$ry,$width,$rect_height,$fill_col,$stroke_col,0);

			($line_start_x,$line_start_y) = ($x,$y+$height+$rect_height);
		}
		elsif( $attr eq "p_end" || $attr eq "q_end" ){
			&svg_ellipse($x+$width/2,$y+$rect_height,$width/2,$height,$fill_col,$stroke_col,1,1);
			&svg_line($line_start_x,$line_start_y,$x,$y+$rect_height,$stroke_col,$stroke_width);
			&svg_line($line_start_x+$width,$line_start_y,$x+$width,$y+$rect_height,$stroke_col,$stroke_width);
			&svg_rect( $x,$y,$rx,$ry,$width,$rect_height,$fill_col,$stroke_col,0);
		}
	}
	foreach my $pos1(sort {$a <=> $b} keys(%{$h_cyto{$chr}})){
		my ($x,$y,$width,$height) = &get_chr_pos($chr,$pos1,$h_cyto{$chr}{$pos1}{"end"});
		my $fill_col = $h_cytoband_color{ $h_cyto{$chr}{$pos1}{"color"} };
		if( $h_chr_attr{$chr}{"p_start"} == $pos1 || $h_chr_attr{$chr}{"q_start"} == $pos1
		 ||  $h_chr_attr{$chr}{"p_end"} == $pos1 || $h_chr_attr{$chr}{"q_end"} == $pos1 ){
		}
		else{
			&svg_rect( $x,$y,$rx,$ry,$width,$height,$fill_col,$stroke_col,0);
			@arr_last = ( $x,$y,$rx,$ry,$width,$height,$fill_col,$stroke_col,0);
			#&svg_rect( $x,$last_y,$rx,$ry,$width,$height,$fill_col,"white",0);
		}

	}
}
