#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2009-05-06
#Modified:
#Description: É¾³ýmatched reads 
my $version=1.00;

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"i=s","o=s","h");
if (!(defined $opts{i} and defined $opts{o} ) || defined $opts{h}) { #necessary arguments
&usage;
}

my $filein=$opts{'i'};
my $fileout=$opts{'o'};
#edit by JC on 2015/12/25
my $fileout_2=$fileout.".pairs";
my $fileout_3=$fileout_2.".example";
my $count=0;
open OUT_2,">$fileout_2";
open OUT_3,">$fileout_3";

#=========================
open IN,"<$filein"; #input file  
open OUT,">$fileout"; #output file  
while (my $aline=<IN>) {

	if($aline=~/>>/){
		$count++;
		print OUT $aline;
		my @temp=split /\s+/,$aline;
		$temp[0]=~s/^>>//g;
		$temp[1]=~s/:[^:]+$//g;
		print OUT_2 "$temp[0]\t$temp[1]\n";
		print OUT_3 "$temp[0]\t$temp[1]\n" if $count<=15;
		}
}

close IN;
close OUT;
close OUT_2;
close OUT_3;
sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i -o
options:
-i input file
-o output file
-h help
USAGE
exit(1);
}

