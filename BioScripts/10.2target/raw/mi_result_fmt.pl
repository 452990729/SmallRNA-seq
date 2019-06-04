use warnings;
use strict;

my $in_flie=shift;
my $out_file=shift;

$/="Complete";

open IN,"<$in_flie";
open TO,">$out_file";

while(<IN>){
	chomp;
	if(/No Hits Found above Threshold/){
		next;}
	else{
		print TO "$_","Complete";
	}
}
close(IN);
close(TO);
