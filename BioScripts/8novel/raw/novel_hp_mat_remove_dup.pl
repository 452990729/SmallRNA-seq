use warnings;
use strict;
use Bio::SeqIO;

die "perl $0 predict.result.csv predict_hairpin.fa predict_mature.fa predict_hairpin.pos" unless(@ARGV==4);

my $input=shift;
my $out_hp=shift;
my $out_mat=shift;
my $out_pos=shift;

open(IN,"$input");
open(OUT1,">$out_hp");
open(OUT2,">$out_mat");
open(OUT3,">$out_pos");
my $marker=0;
my %seq_marker;
while(<IN>){
	chomp;
	if(/novel miRNAs predicted by miRDeep2/){
		<IN>;
		$marker=1;
	}elsif($marker==0){
		next;
	}else{
		my $id="novel_".$marker++;
		my @tmp=split /\t/;
		if(!defined($seq_marker{$tmp[13]})){
			$seq_marker{$tmp[13]}=1;
			print OUT1 ">$id\n$tmp[15]\n";
			print OUT2 ">$id\n$tmp[13]\n>$id*\n$tmp[14]\n";
			print OUT3 "$id\t$tmp[16]\n";
		}
	}
}
close(OUT1);
close(OUT2);
close(OUT3);
