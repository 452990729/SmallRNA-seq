use warnings;
use strict;

my $input=shift;
my $stafile=shift;  # as output
my $uniqfile=shift; # as output 

my %hash;
#total mapped sta
my ($uniqreads,$uniqbases,$totalreads,$totalbases);

#sense mapped sta
my ($seuniqreads,$seuniqbases,$setotalreads,$setotalbases);

#antisense mapped sta
my ($anuniqreads,$anuniqbases,$antotalreads,$antotalbases);

open(STA,">$stafile");
open(UNIQ,">$uniqfile");

open(IN,$input)||die "can't open $input";
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $id=$tmp[0];
	$id=~/^(\S+_\d+\_x\d+)\_\d+$/;
	$id=$1;
	$hash{$id}++;
	my $seq=$tmp[4];
	$totalreads++;
	$totalbases+=length($seq);

	if($tmp[1] eq "-"){
		$seq = reverse($seq);
		$seq =~ y/CGATUcgatu/GCTAAGCTAA/;
		$antotalreads++;
		$antotalbases+=length($seq);
	}else{
		$setotalreads++;
		$setotalbases+=length($seq);
	}
	if($hash{$id}==1){
		$uniqreads++;
		$uniqbases+=length($seq);;
		print UNIQ ">",$id,"\n",$seq,"\n";
		if($tmp[1] eq "-"){
			$anuniqreads++;
			$anuniqbases+=length($seq);
		}else{
			$seuniqreads++;
			$seuniqbases+=length($seq);
		}	
	}
}
close(IN);


print STA "Statu\tTotal reads\tTotal bases (bp)\tUniq reads\tUniq bases (bp)\n";
print STA "Total Mapped small RNA\t$totalreads\t$totalbases\t$uniqreads\t$uniqbases\n";
print STA "Total Sense Mapped small RNA\t$setotalreads\t$setotalbases\t$seuniqreads\t$seuniqbases\n";
print STA "Total Antisense Mapped small RNA\t$antotalreads\t$antotalbases\t$anuniqreads\t$anuniqbases\n";

