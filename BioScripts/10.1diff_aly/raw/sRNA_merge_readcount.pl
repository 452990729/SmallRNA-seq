use strict;
use warnings;

sub usage{
print STDERR <<USAGE;
=========================================================================
example: perl $0 known/mature.readcount,novel/mature.readcount outdir/out_prefix

readcount file format:
*****************************
type\tsample1\tsample2\tsample3...\n
id\t1\t2\t3...\n
*****************************
=========================================================================
USAGE
exit 0;
}

&usage unless (@ARGV==2);
my @infile=split(",",shift);
my $out_prefix=shift;

my (@sample,$type);
my %count_type_sample;
my %total_count;
foreach my $i (@infile){
	open(IN,"$i");
	my $header=<IN>;
	chomp($header);
	my @tmp_header=split("\t",$header);
	if( !($type) ) {
		($type,@sample)=split("\t",$header);
	}
	while(<IN>){
		chomp;
		my @tmp=split /\t/ ;
		foreach my $j(1..$#tmp){
			$count_type_sample{$tmp[0]}{$tmp_header[$j]}=$tmp[$j];
			$total_count{$tmp_header[$j]}+=$tmp[$j];
		}
	}
	close(IN);
}

open(OUT1,">$out_prefix.readcount");
open(OUT2,">$out_prefix.tpm");
print OUT1 "sRNA\t",join("\t",@sample),"\n";
print OUT2 "sRNA\t",join("\t",@sample),"\n";
#print $type,"\t",join("\t",@sample),"\n";
foreach my $i(sort {$a cmp $b} keys %count_type_sample){
	print OUT1 $i;
	print OUT2 $i;
	foreach my $j(@sample){
		print OUT1 "\t",$count_type_sample{$i}{$j};
		print OUT2 "\t",(1000000*$count_type_sample{$i}{$j})/$total_count{$j};
	}
	print OUT1 "\n";
	print OUT2 "\n";
}

