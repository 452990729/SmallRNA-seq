#!/usr/bin/perl -w
use strict;
use Getopt::Long;
sub usage
{
	print STDERR <<USAGE;
=============================================================================================
Description	format split blast out 

Version:1.0
2013-07-30 zhangyingying\@novogene.cn

Options
<Requried>
	-xml <s>	:	blast result: *.xml
	-geneid <s>	:	diff gene ID
	-out <s>	:	outfile: diffgeneID.blastout
=============================================================================================
USAGE
}
my($help,$xml,$geneid,$out);
GetOptions(
	"h|?|help"=>\$help,
	"geneid=s"=>\$geneid,
	"xml=s"=>\$xml,
	"out=s"=>\$out,
);

if(defined($help) || !defined($xml) || !defined($geneid) || !defined($out)){
	&usage;
	exit 0;
}
##########################################
#get diff gene id
my %hash;
open GENE, "$geneid";
while(<GENE>){
	chomp;
	$hash{$_}=1;
}
close GENE;
##########################################
#split xml
open XML,"$xml";
open OUT,">$out";
my $tag=0;
$/ = "<Iteration>";
my $header=<XML>;
chomp($header);
print OUT $header;#header
while (my $line=<XML>){
	chomp($line);
	if($line=~/Iteration_query-def\>(.*?)\<\/Iteration_query-def/){
		if(defined($hash{$1})){
			if($line=~ /BlastOutput_iterations/){
				$tag=1;
			}
			print OUT "<Iteration>$line";
		}
	}
}
$/ = "\n";
if($tag eq "0"){
	print OUT "</BlastOutput_iterations>\n</BlastOutput>\n";
}
close XML;
close OUT;

