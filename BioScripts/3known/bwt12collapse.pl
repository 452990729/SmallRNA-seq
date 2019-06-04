use warnings;
use strict;

my $input=shift;

open(IN,$input)||die "can't open $input";
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $id=$tmp[0];
	my $seq=$tmp[4];
	if($tmp[1] eq "-"){
		$seq = reverse($seq);
		$seq =~ y/CGATUcgatu/GCTAAGCTAA/;
	}
	print ">",$id,"\n",$seq,"\n";
}
close(IN);

