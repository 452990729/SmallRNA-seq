use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;

my $usage=<<END;
Program:
	this program sumarize novel miRNA
History:
	2015/5/27	writer jc First release
	2015/10/28	edit by jc

Eg:
perl $0 -readcount	mature.readcount 
	-h	hairpin.fa 
	-m	mature.fa 
	-str	hairpin.str (predict hirpin second structure)
	-pos	hairpin.pos (predict hirpin position)
	-o	novel_hairpin_mature.xls (output file)

output: novel_hairpin_mature.xls
END

my ($csv,$hirpin,$mature,$out,$hirpin_str,$hirpin_pos);

GetOptions(
	"readcount=s" =>\$csv,
	"h=s"	=>\$hirpin,
	"m=s"	=>\$mature,
	"str=s"	=>\$hirpin_str,
	"pos=s" =>\$hirpin_pos,
	"o=s"	=>\$out,
);

if ((! defined $csv) or (! defined $hirpin) or ( ! defined $mature)or (!defined $hirpin_str)){
	die $usage;
}

if (defined $out){
	open OUT,">$out";
}else{
	open OUT,">novel_hairpin_mature.xls";
}

open CSV,"$csv" or die "can't open $csv";
my $in_hirpin=new Bio::SeqIO(-file =>$hirpin);
my $in_mature=new Bio::SeqIO(-file =>$mature);

#get hairpin_mature.xls's head==========================================

my $head="precursor_id\tprecursor_seq\tmfe\tposition\tmature_id\tmature_seq\t";
chomp(my $csv_head=<CSV>);
my @csv_head=split /\t/,$csv_head;
$head .=(join "\t",@csv_head[1..$#csv_head])."\n";	#get sample id

print OUT "$head";
#======================================================================

#get hirpin/mature id sequence ========================================
my %hash_hirpin;	#input hirpin.fa ;
my %hash_mature;	#input mature.fa ;

while (my $seq=$in_hirpin -> next_seq()){
	$hash_hirpin{$seq ->id}=$seq -> seq;
}
while (my $seq=$in_mature -> next_seq()){
	$hash_mature{$seq -> id}=$seq -> seq;
}

#=======================================================================
##stat mfe of predict hirpin ===============================

open STR,"$hirpin_str";
my %hash_str;
while (<STR>){
        chomp;
        if (/^>(\S*)/){
                my $hirpin_id=$1;
                my $hirpin_seq=<STR>;
                chomp (my $hirpin_stucture=<STR>);
                my ($mfe)= $hirpin_stucture =~ / \((.+?)\)$/;
                $hash_str{$hirpin_id}=$mfe;
        }
}
close STR;
#============================================================

#stat hirpin position =======================================
open POS,"$hirpin_pos";
my %hash_pos;
while (<POS>){
	chomp;
	$hash_pos{(split /\t/)[0]}=(split /\t/)[1];
}
close POS; 
#============================================================

#out put=====================================================


while (<CSV>){
	chomp;
	my @temp=split /\t/;
		my $sample_readcount=(join "\t",@temp[1..$#temp]);
		print OUT "$temp[0]\t$hash_hirpin{$temp[0]}\t$hash_str{$temp[0]}\t$hash_pos{$temp[0]}\t$temp[0]\t$hash_mature{$temp[0]}\t$sample_readcount\n";
			
}
close CSV;
close OUT;
#====================================================================

