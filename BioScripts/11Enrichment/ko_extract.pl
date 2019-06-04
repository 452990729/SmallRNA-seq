#!/usr/bin/perl
die "perl $0 <gene list> <ko annotation> <out>\n" unless @ARGV==3;

$list=shift;
$ko=shift;
$out=shift;
my %ko_hs=();

open FILE,$ko;

<FILE>;   ##skipp first line,  and feed in Species(ko) to be compatable with new version kobas2.0-20140801
my $head = "##ko\n";
$head.=<FILE>;
$head.=<FILE>;
$head.=<FILE>;

while(<FILE>){
	chomp;
	my @array=split/\t/;
	$ko_hs{$array[0]}=$array[1];

	if(/--------------------/){
		last;
	}
}
close FILE;

open OUT,">$out";
print OUT $head;
open LIST,$list;
while(<LIST>){
	chomp;
	print OUT $_."\t".$ko_hs{$_}."\n";
}
close OUT; close LIST;
