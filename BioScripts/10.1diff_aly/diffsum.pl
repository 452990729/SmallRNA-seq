#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

sub usage
{
        print STDERR <<USAGE;
==================================================================================
Description  check diff numbers
Options
<Optional>
        -diffdir    <s>  :   the diff aly dir, default pwd/Diff_TR
	-cutoff     <i>  :   the cutoff value for diff number,default 10
===============================================================================
USAGE
}
my($diff,$cutoff,$help);
GetOptions(
        "h|?|help"=>\$help,
                "diffdir=s"  =>\$diff,
                "cutoff=i"   =>\$cutoff,
);

if(defined($help)){
        &usage;
        exit 0;
}

$cutoff||=10;
my $pwd = `pwd`;
$pwd =~ s/\s+//g;
$diff ||= "$pwd/diffAnalysisResult";
my $out = "$diff/diff_sum.txt";
my @groups = `ls $diff/ | grep vs | less`;
my $group;
my @diffnums;
push(@diffnums,"GROUP\tDIFF\tUP\tDOWN\n");
foreach $group(@groups)
{
    $group =~ s/\s+//g;
    if ($group!~ m/vs.*vs/)
    {
	my $diff_count = `wc -l $diff/$group/*DElist.txt`;
	my $up_count = `wc -l $diff/$group/*DElist_up.txt`;  
	my $down_count = `wc -l $diff/$group/*DElist_down.txt`;  
	my @diff_count_split = split /\s+/,$diff_count;
	my @up_count_split = split /\s+/,$up_count;
	my @down_count_split = split /\s+/,$down_count;
	push(@diffnums,"$group\t$diff_count_split[0]\t$up_count_split[0]\t$down_count_split[0]\n");
	if($diff_count_split[0]<$cutoff||$up_count_split[0]<$cutoff||$down_count_split[0]<$cutoff){
		print "###################!!!!!WARNINGS!!!!!#########################\n$group\'s diffgene numbers are too low!(<$cutoff)\n";
		print "$diff_count_split[0]\t$up_count_split[0]\t$down_count_split[0]\n"
	}
    }
}
open(OUT,">$out");
print OUT @diffnums;
close OUT;
