use strict;
use warnings;

open IN,shift;
my $id;
my %hash;
my %hash_sample;
while(<IN>){
	chomp;
	if (/^>(\S+)$/){
		$id=$1;
		$hash{$id}=1;
	}else{
		if (/^(\S+)_\d+_x\d+/){
			$hash_sample{$1}{$id}=1;
		}
	}
}

print "Total\t",join("\t",sort keys %hash_sample),"\n";
my $all=keys %hash;
print "$all\t";
for (sort keys	%hash_sample){
	my $num=keys $hash_sample{$_};
	print "$num\t";
}
print "\n";
close IN;

__DATA__
open IN_1,shift;
open IN_2,shift;

my (%hash,%hash_sample,%hash_pre);
while (<IN_1>){
	chomp;
	$hash{(split /\t/)[2]}=1;
}

close IN_1;

while (<IN_2>){
	chomp;
	my @temp=split /\t/;
	my ($sample)=$temp[0]=~ /^(\S\S\S)/;
	next unless exists $hash{$temp[2]};
	$hash_pre{$temp[2]}=1;
	$hash_sample{$sample}{$temp[2]}+=1;
}

my @all_pre=keys %hash_pre;
print "all_pre\t",$#all_pre+1,"\n";
for my $sample(sort keys %hash_sample){
	my @each_sample_pre=keys $hash_sample{$sample};
	print "$sample\t",$#each_sample_pre+1,"\n";
}
close IN_2;
