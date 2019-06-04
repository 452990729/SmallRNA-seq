#!/usr/bin/perl -w
use warnings;
use strict;

if (@ARGV < 3)
{
	print "Usage: $0 directoryGeneSeq selectInfo directoryOutput\n";
	exit;
}
my ($dirGeneSeq, $dirInfo, $dirOutput) = @ARGV;
my %diffGene;
my @geneIDs;
#$dirOutput ||= "/database/rna/software/KOBAS/output/20120419/test_result_2.fa";

open INPUT, "<$dirGeneSeq";
open INFO, "<$dirInfo";
open OUTPUT, ">$dirOutput";

while (<INFO>)
{
	chomp;
	my @tmp=split;
	$diffGene{$tmp[0]} = 1;
}
close INFO;
#print "$geneIDs[23]";

my $id;
while (<INPUT>)
{
	if (/^>(\S+)/)
	{
		$id=$1;
		if ($diffGene{$1})
		{
			$diffGene{$1} = <INPUT>;
			print OUTPUT ">$1\n$diffGene{$1}";
		}
	}elsif($diffGene{$id}){
	print OUTPUT $_;
	}
}
close INPUT;
close OUTPUT;
