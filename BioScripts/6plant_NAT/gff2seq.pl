#!/bin/perl -w
use strict;
use Getopt::Long;
my ($help, $gff, $feature, $seq, $out);
GetOptions(
	"h|help"		=>\$help,
	"g|gff=s"		=>\$gff,
	"f|feature=s"	=>\$feature,
	"s|seq=s"		=>\$seq,
	"o|outfile=s"	=>\$out
);
my $usage=<<END;
-----------------------------------------------------
perl $0 -g *.gff -f feature -s *.fa -o outfile

-h|help		help
-g|gff=s	gff file
-f|feature=s 	3rd column eg. gene/exon ...
-s|seq=s	sequence file
-o|outfile=s	output file
-----------------------------------------------------
END
die $usage if ($help or !$gff or !$seq);

open S, $seq or die "Can't open $seq: $!";
my(%hash, $id);
while(<S>)
{
	chomp;
	if(/^>(\S+)/)
	{
		$id=$1;
	}
	else
	{
		$hash{$id}.=$_;
	}
}
close S;

open G, $gff or die "Can't open $gff: $!"; 	
#Chr1	TAIR10	chromosome	1	30427671	.	.	.	ID=Chr1;Name=Chr1
#Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
#Chr1	TAIR10	mRNA	3631	5899	.	+	.	ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
#Chr1	TAIR10	protein	3760	5630	.	+	.	ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
open O, ">$out" or die "Can't output $out: $!";
my $fa;
while(<G>)
{
	chomp;
	my ($chr, $source, $type, $start, $end, $score, $strand, $frame, $attr)=split;
	#my $id=$1 if ($attr=~/ID=(\w+);+?/);
	my $id=$1 if ($attr=~/ID=(\S+?);/);
	if($type eq $feature)
	{
		if($strand eq "+" && defined $hash{$chr})
		{
			$fa=substr($hash{$chr}, $start-1, $end-$start+1);
			print O ">$id\n$fa\n";
		}
		elsif($strand eq "-" && defined $hash{$chr})
		{
			$fa=substr($hash{$chr}, $start-1, $end-$start+1);
			$fa=~tr/ATCGatcg/TAGCtagc/;
			$fa=reverse $fa;
			print O ">$id\n$fa\n";
		}
	}
}
close G;
 
			
