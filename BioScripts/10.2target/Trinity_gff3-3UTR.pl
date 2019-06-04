use warnings;
use strict;

#usge: perl region2UTR.pl $refer $gff3file >3'UTR_file
use Bio::Perl;
use Bio::DB::Fasta;
use List::Util qw(reduce);

my $db;
my $refer = shift;
my $cds_gff=shift;
my $min_len=shift;
if(-e "$refer.index"){
	`rm -rf $refer.index`;
}
$db = Bio::DB::Fasta->new("$refer");

my $seq=();
open(FILE,$cds_gff);
while(<FILE>){
	chomp;
if(/^\S/){
	my @array=split /\t/;
	if($array[2] eq "CDS"){
		$seq=$db->subseq($array[0]);
		chomp($seq);
		my $end=length($seq);
		if($array[6] eq '+'){
			$end=$end;
			if($end-$array[4]>=$min_len){
				$seq=$db->subseq($array[0],$array[4] + 1,$end);
				print ">","$array[0]:",$array[4] + 1,"-",$end,"(+) type:3-UTR len:",$end-$array[4],"\n",$seq,"\n";
#				print ">",$array[0]," type:3-UTR len:",$end-$array[4]," (+) ","$array[0]:",$array[4] + 2,"-",$end+1,"(+)\n",$seq,"\n";
			}
		}
		if($array[6] eq '-'){
			if($array[3]>$min_len){
				$seq=$db->subseq($array[0],1,$array[3] - 1);
				$seq=~tr/ATCG/TAGC/;
				$seq=reverse $seq;
				print ">$array[0]:1-",$array[3]-1,"(-) type:3-UTR len:",$array[3]-1,"\n",$seq."\n";
#				print ">",$array[0]," type:3-UTR len:$array[3] (-) $array[0]:1-$array[3](-)\n",$seq."\n";
			}
		}
	}
}
}
close(FILE);


