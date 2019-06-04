#By zhangyu, 2013.01.10
use strict;
use warnings;

die "perl $0 <in_pairs> <inAnnot>" unless (@ARGV==2);

my $in=shift;
my $anno=shift;
my $out=$in.".mg";
my $out2=$in.".annotate";

#if($anno =~ /Blast_NT.xls/){
#	my $anno1="Blast_NT.hit1";
#	`cut -f1,16 $anno >$anno1`;
#	$anno=$anno1;
#}

open IN,"<$in";
open IN2,"<$anno";
open OUT,">$out";

my %hash;
my %hash2;

while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $geneID=$tmp[1];
	my $miRNA=$tmp[0];
	${$hash{$geneID}{$miRNA}}=0;
}
foreach my $key(keys %hash){
	print OUT "$key\t";
	#foreach my $key2(keys %{$hash{$key}}){
	#	print OUT "$key2,";
	#}
	print OUT join(",", keys %{$hash{$key}});  # 2014-2-22 wangshaobin
	print OUT "\n";
}
close IN;
close OUT;
open OUT, "<$out";
while(<OUT>){
        chomp;
        my @tmp2=split /\t/;
        my $geneid=$tmp2[0];
		my $miRNAid=$tmp2[1];
        $hash2{$geneid}=$miRNAid;
}
close OUT;

#open OUT2,">>$out2";
open OUT2,">$out2";
my $headline=<IN2>;
chomp($headline);
my @header=split /\t/,$headline;
shift @header;
my $head=join("\t",@header);
print OUT2 "geneID\tmiRNA\t$head\n";
 
while(<IN2>){
	chomp;
	my @atm=split /\t/;
	my $id=shift @atm;
	my $anno=join("\t",@atm);
	if(exists $hash2{$id}){
		print OUT2 "$id\t$hash2{$id}\t$anno\n";
	}
}
close IN2;
close OUT2;

