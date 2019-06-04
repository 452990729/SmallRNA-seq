use warnings;
use strict;
use Bio::SeqIO;

die "perl $0 hairpin.fa mature.fa hairpin_out.fa mature_out.fa" unless(@ARGV==4);

my $in_hp=shift;
my $in_mat=shift;
my $out_hp=shift;
my $out_mat=shift;

open(OUT1,">$out_hp");
open(OUT2,">$out_mat");
my %seq_mat;
my %hp_marker;

#read in_hp file
my %hp_seq;
my $in1 = new Bio::SeqIO(-file => $in_hp);
while(my $seq = $in1->next_seq()){
	my $sequence=$seq->seq;
	$sequence=~s/[^AGCTUCNagctun]/N/g;
	$hp_seq{$seq->id}=$sequence;
}

#read in_mat file
my $in2 = new Bio::SeqIO(-file => $in_mat);
while(my $seq = $in2->next_seq()){
	my $sequence=$seq->seq;
	$sequence=~s/[^AGCTUCNagctun]/N/g;
	if(!(defined($seq_mat{$sequence})||defined($seq_mat{substr($sequence,1)})||defined($seq_mat{substr($sequence,2)})||defined($seq_mat{substr($sequence,1,length($sequence)-2)})||defined($seq_mat{substr($sequence,0,length($sequence)-2)})||defined($seq_mat{substr($sequence,0,length($sequence)-1)}))){
		$seq_mat{$seq->seq}=$seq->id;
		my $id=$seq->id;
		$id =~ s/\-5p//g;
		$id =~ s/\-3p//g;
		$id =~ s/\*//g;
		$id =~ s/\.\d+//g;
		foreach my $id1(keys %hp_seq){
			if(($id1 =~ /$id/i) || ($id =~ /$id1/i)){
				if(!defined($hp_marker{$id1})){
					print OUT1 ">",$id1,"\n",$hp_seq{$id1},"\n";
					$hp_marker{$id1}=1;
				}
			}
		}
		print OUT2 ">",$seq->id,"\n",$sequence,"\n";
	}
}
