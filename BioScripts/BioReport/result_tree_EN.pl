#!usr/bin/perl -w
use strict;
use FindBin '$Bin';
open(IN,"$Bin/sRNA_tree_template_EN.xls");
my %hash;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	$hash{$tmp[0]}=$tmp[1];
}
close(IN);
open(IN1,"$ARGV[0]") or die "Usage: perl result_tree.pl DirectoryTree.html DirectoryTree.TL.html";
open(OUT,">$ARGV[1]")or die "Usage: perl result_tree.pl DirectoryTree.html DirectoryTree.TL.html";
$/="<\/a>";
while(<IN1>){
	chomp;
	if(/<br>/){
		if(/(.*\d\.)(.*)$/){
			#print "$1\n$2\n";die;
			if(defined $hash{$2}){
				print OUT "$1\t$2:\t$hash{$2}\t<\/a>\n";
			}else{
				print OUT "$_<\/a>\n";
			}
		}elsif(/<br><br>/){
			print OUT "$_\n";
		}
	}else{
		print OUT "$_<\/a>\n";
	}
}
close IN1;
close OUT;
`sed -i s/charset=iso-8859-1/charset=UTF-8/g $ARGV[1]`;
