use warnings;
use strict;

die "usage: check 3'utr or gene.fa,input gtf and 3'utr(or gene.fa)" unless(@ARGV ==2);

open GTF,"<$ARGV[0]";
open UTR,"<$ARGV[1]";

my %hash;
while (<GTF>){
	chomp;
	next if (/^#/);
	my @temp=split /\t/;
	if ($temp[8]=~ /transcript_id \"(\S+?)\";/){
	$hash{$1} +=1;
	}
}

while (<UTR>){
	chomp;
	s/\s*$//g;
	next if (! /^>/);
	if (/^>(.*)/){
		if (not exists $hash{$1}){
			die "the 3'utr/gene.fa is not ok in the 3'utr/gene.fa's $1\n";
			exit (1);
		}
	}
}
print "The 3'utr or gene.fa is ok\n";
close GTF;
close UTR;
