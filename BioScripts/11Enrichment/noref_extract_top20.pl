#!/usr/bin/perl
use POSIX qw(log10);
my $file=shift;
open IN, $file, or die;
<IN>;
<IN>;
<IN>;
<IN>;
<IN>;

#print "pathway_term\tRF\tQvalue\t-log10\(Qvalue\)\n";
print "pathway_term\trich_factor\tqvalue\tgene_number\n";

my $i=0;
while (<IN>) {
	chomp;
    next if(/^@@/);
	$i++;
	if ($i>20 || /--------/ || /^$/) 
	{
		last;
	}
	else
	{
		chomp;
		my @tmp=split /\t/;
		my $RF=$tmp[3]/$tmp[4];
		#my $Qvalue=-log10($tmp[6]);
		$tmp[0]=~s/'/"'"/;
		print "$tmp[0]\t$RF\t$tmp[6]\t$tmp[3]\n";
	}
}
close IN;
