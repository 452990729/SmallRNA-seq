#Edit by jiangxiaoxue
###2012/09/25

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Config::Tiny;

sub usage{
	print STDERR <<USAGE;
======================================================================================================================================
Usage: perl $0 [options]
options:
	-h|?|--help	help information
	-r rfam.fasta	Mapping reference in fasta format(bwtindex prefix=fastaname)
			[default="/PUBLIC/database/Common/Rfam/dna/Rfam.dna.fasta"]
	-i read.fa	not special format
	-s sampid	sample name
	-o outdir       output project
	-T [int]	number of alignment threads to launch during bowtie mapping(default:10)
======================================================================================================================================
USAGE
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $bowtie = $Config->{software}->{bowtie1};
my $Rfam_2017 = $Config->{database}->{Rfam_2017}; 

my ($help,$refer,$read,$sampid,$outdir,$threads);

GetOptions(
	"h|?|help"=>\$help,
	"r=s"=>\$refer,
	"i=s"=>\$read,
	"s=s"=>\$sampid,
	"o=s"=>\$outdir,
	"T=i"=>\$threads
);

$threads ||= 10;
$refer ||= $Rfam_2017;
if(!defined($read) || !defined($sampid)  || defined($help)){
	&usage;
	exit 0;
}

$outdir ||= `pwd`;
$outdir =~ s/\s+//g;


open(REF,"$refer");
my %gi_spe;
while(<REF>){
	if(/^>/){
		/^>(\S+)/;
		my $gi=$1;
		/^>.*:(.*)/;
		my $spe=$1;
		$gi_spe{$gi}=$spe;
	}
}


my %count;

`$bowtie/bowtie -v 0 -p $threads -f $refer $read >$outdir/$sampid\_category.bwt 2>$outdir/$sampid\_category.bwtlog`;
my $total=`grep -c ">" $read`;
chomp($total);
my $map;
my %type_count;
open(BWT,"$outdir/$sampid\_category.bwt");
open(OUT1,">$outdir/$sampid\_category_detail.txt");
while(<BWT>){
	chomp;
	my @tmp=split /\t/;
	$map++;
	print OUT1 "$tmp[0]\t$tmp[4]\t$tmp[1]\t",$gi_spe{$tmp[2]},"\t1\n";
	$count{$gi_spe{$tmp[2]}}++;
	my @tmp1=split(";",$tmp[2]);
	$type_count{$tmp1[1]}++;
}
close(BWT);
close(OUT1);

open(OUT2,">$outdir/$sampid\_category.types");
print OUT2 "Rfam_type\tmapped_sRNA\tmapped_sRNA/total_sRNA(%)\tmapped_sRNA/total_mapped_sRNA(%)\n";
foreach my $i(sort {$type_count{$b}<=>$type_count{$a}} keys %type_count){
	my $rate1=sprintf("%.2f",$type_count{$i}*100/$total);
	my $rate2=sprintf("%.2f",$type_count{$i}*100/$map);
	print OUT2 $i,"\t",$type_count{$i},"\t",$rate1,"\t",$rate2,"\n";
}
close(OUT2);

open(OUT,">$outdir/$sampid\_category.log");
print OUT "total_sRNA=>\t",$total,"\t100\n";
print OUT "mapped_sRNA=>\t",$map,"\t100\n","mapped_sRNA/total_sRNA(%)=>\t",sprintf("%.2f",$map*100/$total),"\n\n";
print OUT "type\tspe_mapped_sRNA\tpercent_total(%)\tpercent_mappes(%)\n";
foreach my $key ( sort { $count{$b} <=> $count{$a} } keys %count ) {
	my $rate1=sprintf("%.2f",$count{$key}*100/$total);
	my $rate2=sprintf("%.2f",$count{$key}*100/$map);
	print OUT $key,"\t",$count{$key},"\t",$rate1,"\t",$rate2,"\n";
}
close(OUT);

