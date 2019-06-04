use warnings;
use strict;

die "perl $0 prefix_file miRBase.fa" unless (@ARGV == 2);

my %hash=();
my $pre_file=shift;
my $infile=shift;

open(PRE,"$pre_file")||die "can't open $pre_file\n";
while(<PRE>){
	chomp;
	$hash{$_}=1;
}
close(PRE);

my $marker=0;
open(FILE,$infile);
while(<FILE>){
	chomp;
	my $tmp=$_;
	if(/^>/){
		s/\-.*//;
		s/^>//;
		if(defined($hash{$_})){
			$marker=1;
			print "$tmp\n";
		}
		else{
			$marker=0;
		}
	}
	elsif($marker){
		print "$tmp\n";
	}

}
close(FILE);
