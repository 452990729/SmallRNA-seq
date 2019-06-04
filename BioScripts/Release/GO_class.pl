#!/usr/bin/perl -w
use strict;

die "Usage:perl $0 <in.GO> <in.gene_ontology.1_2.obo.ab> <out.GO.class>" unless(@ARGV==3);

my $in=shift;
my $ab=shift;
my $out=shift;

open IN,"<$in";
open AB,"<$ab";
open OUT,">$out";
my @goterm=qw/biological_process cellular_component molecular_function/;

my %GO;
while(<AB>){
	chomp;
	my @line=split /\t/,$_;
	@{$GO{$line[0]}}=($line[2],$line[1]);
}
my %gene;
my %class;
my %des;
while(<IN>){
	chomp;
	my @array=split /\t/,$_;
	foreach(my $i=1;$i<=$#array;$i++){
		if(defined(${$GO{$array[$i]}}[0])){
			push @{$gene{$array[0]}{${$GO{$array[$i]}}[0]}},$array[$i];
		}else{
			print "no description:\t$array[0]\t$array[$i]\n";
		}
	}
	foreach(my $j=0;$j<=$#goterm;$j++){
		if(defined(${$gene{$array[0]}{$goterm[$j]}}[0])){
			$class{$array[0]}{$goterm[$j]}=join("//",@{$gene{$array[0]}{$goterm[$j]}});
			my @tmp;
			foreach(my $k=0;$k<=$#{$gene{$array[0]}{$goterm[$j]}};$k++){
				$tmp[$k]=${$GO{${$gene{$array[0]}{$goterm[$j]}}[$k]}}[1];
			}
			$des{$array[0]}{$goterm[$j]}=join("//",@tmp);
		}else{
			$class{$array[0]}{$goterm[$j]}="--";
			$des{$array[0]}{$goterm[$j]}="--";
		}
	}	
}

print OUT "Gene_ID\t$goterm[0]\t$goterm[0]\_description\t$goterm[1]\t$goterm[1]\_description\t$goterm[2]\t$goterm[2]\_description\n";
foreach my $key1(sort keys %gene){
	print OUT "$key1\t$class{$key1}{$goterm[0]}\t$des{$key1}{$goterm[0]}\t$class{$key1}{$goterm[1]}\t$des{$key1}{$goterm[1]}\t$class{$key1}{$goterm[2]}\t$des{$key1}{$goterm[2]}\n";
}

close IN;
close AB;
close OUT;
