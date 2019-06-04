#By zhangyu, 2013.01.10
use strict;
use warnings;

die "perl $0 <in_pairs> <inAnnot>" unless (@ARGV==2);

my $in=shift;
my $anno=shift;
my $out=$in.".mg";

open IN,"<$in";
open IN2,"<$anno";
open OUT,">$out";

my %hash;

my $headline=<IN2>;
chomp($headline);
my @header=split /\t/,$headline;
shift @header;
my $head=join("\t",@header);
print OUT "miRNA\tgeneID\t$head\n";
 
while(<IN2>){
	chomp;
	my @atm=split /\t/;
	my $id=shift @atm;
	my $anno=join("\t",@atm);
	$hash{$id}=$anno;
}

while(<IN>){
	chomp;
	my @tmp=split / /;
	print OUT "$_\t$hash{$tmp[1]}\n";
}
