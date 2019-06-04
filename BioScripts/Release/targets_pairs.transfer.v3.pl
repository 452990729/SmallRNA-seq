#By zhangyu, 2013.01.10
use strict;
use warnings;

die "perl $0 <in> <out1>" unless (@ARGV==2);

my $in=shift;
my $out1=shift;

open IN,"<$in";
open OUT1,">>$out1";

my %hash1;
my %hash2;

while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $miRNA=$tmp[0];
	my $mRNA=$tmp[1];
	$hash1{$miRNA}{$mRNA}=0;
	$hash2{$mRNA}{$miRNA}=0;
}
close IN;

print OUT1 "ID\tTargetCorrespond\n";

foreach my $key(keys %hash1){
	print OUT1 "$key\t";
	foreach my $key2(keys %{$hash1{$key}}){
		print OUT1 "$key2,";
	}
	print OUT1 "\n";
}
foreach my $key3(keys %hash2){
        print OUT1 "$key3\t";
        foreach my $key4(keys %{$hash2{$key3}}){
                print OUT1 "$key4,";
        }
        print OUT1 "\n";
}
close OUT1;

