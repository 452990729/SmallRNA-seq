use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my $usage=<<END;
Program:
	This program get hairpin and mature information from mirDeep quantifier result miRNAs_expressed_all_samples_xxx.known.csv

History:
	2015/5/27	writer jc First release

Eg:
perl $0 -readcount	mature.readcount 
	-h	hairpin.fa
	-m	mature.fa
	-pairs	hairpin_mature.pairs
	-o	hairpin_mature.xls (output file)

output: known_hairpin_mature.xls
END

my ($csv,$hirpin,$mature,$out,$pairs);

GetOptions(
	"readcount=s" =>\$csv,
	"h=s"	=>\$hirpin,
	"m=s"	=>\$mature,
	"pairs=s" =>\$pairs,
	"o=s"	=>\$out,
);

if ((! defined $csv) or (! defined $hirpin) or ( ! defined $mature)){
	die $usage;
}

if (defined $out){
	open OUT,">$out";
}else{
	open OUT,">known_hairpin_mature.xls";
}

open CSV,"$csv" or die "can't open $csv";
my $in_hirpin=new Bio::SeqIO(-file =>$hirpin);
my $in_mature=new Bio::SeqIO(-file =>$mature);

my $head="precursor_id\tprecursor_seq\tmfe\tposition\tmature_id\tmature_seq\t";
chomp(my $csv_head=<CSV>);
my @csv_head=split /\t/,$csv_head;
$head .=join("\t",@csv_head[1..$#csv_head]);	#get sample id
#===test 
print OUT "$head\n";

my %hash_hirpin;	#input hirpin.fa ;
my %hash_mature;	#input mature.fa ;

while (my $seq=$in_hirpin -> next_seq()){
	$hash_hirpin{$seq ->id}=$seq -> seq;
}
while (my $seq=$in_mature -> next_seq()){
	$hash_mature{$seq -> id}=$seq -> seq;
}

open PAIR,$pairs;
my %pairs;
while(<PAIR>){
	chomp;
	$pairs{(split/\t/)[0]}=(split/\t/)[1];
}
close PAIR;

while (<CSV>){
	chomp;
	my @temp=split /\t/;
	my $result='';
	if ( exists $hash_hirpin{$pairs{$temp[0]}}){
		$result="$pairs{$temp[0]}\t$hash_hirpin{$pairs{$temp[0]}}\t";
	}else{
		die "wrong with $temp[0], not hairpin id in ($!)\n";
	}
	$result .="--\t--\t$temp[0]\t$hash_mature{$temp[0]}\t".join("\t",@temp[1..$#temp])."\n";
	print OUT $result;
}
	
close OUT;			
			
close CSV;

