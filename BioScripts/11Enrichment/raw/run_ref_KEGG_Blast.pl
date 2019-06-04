#Edit by jiangxiaoxue
##2012/11/30
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
	-i <s>	different-miRNA-list file
	-t <s>	miRNA-target-pairs file
	--head	different-gene-list with headinfo
	-species <s>	abbr for species(from /PUBLIC/software/public/Annotation/kobas2.0-20120208/scripts/sp.list)
	-genome	<s>	genome sequece
	-gtf <s>	Gene sets for your organism, must be a gtf file.
	-o <s>	project name(suggest using the different express analysis each group name like sC1_vs_sC2)
	-h|?|help	show this help
======================================================================================================================================
USAGE
}
my ($help,$input,$target,$gtf,$species,$goann,$genome,$head,$com);
GetOptions("h|?|help" => \$help,
	"i=s"=>\$input,
	"t=s"=>\$target,
	"gtf=s"=>\$gtf,
	"species=s"=>\$species,
	"genome=s"=>\$genome,
	"head"=>\$head,
	"o=s"=>\$com,
);

if(defined($help)||!defined($input)||!defined($target)||!defined($gtf)||!defined($species)||!defined($genome)){
	&usage;
	exit(0);
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $genome_path=get_ful_path($genome);
my $gtf_path=get_ful_path($gtf);

$genome=basename($genome);
$gtf=basename($gtf);

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
open(OUT_KEGG,">$pwd/KEGG_Blast.sh");

if(defined($head)){
        `awk '{if(NR>1){print}}' $input |sort -k 1,1 >$pwd/diffmiRNAID`;
}else{
        `sort -k 1,1 $input >$pwd/diffmiRNAID`;
}
`sort -k 1,1 $target|join $pwd/diffmiRNAID - |sort -u >$pwd/diffmiRNA-gene.pairs`;
`awk '{print \$2}' $pwd/diffmiRNA-gene.pairs|sort -u >$pwd/diffmiRNA-geneid`;
`ln -sf $genome_path $pwd/$genome`;
`ln -sf $gtf_path $pwd/$gtf`;
#print OUT_KEGG "echo ================== Run different express miRNA target gene enrichment analysis==================\n";
#print OUT_KEGG "echo Start Time:\ndate\n";

print OUT_KEGG "echo ================== Run KOBAS==================\n";
print OUT_KEGG "echo Start Time:\ndate\n";
print OUT_KEGG "$perlExec $Bin/runKEGG.pl  -directory-gtf $gtf -gene-id diffmiRNA-geneid -exist-genome y -directory-kobas-output $pwd -directory-genome $genome -species $species -samplename $com\n";
print OUT_KEGG "echo End Time:\ndate\n";
print OUT_KEGG "cd $pwd\n";
#close(OUT);
close(OUT_KEGG);
sub get_ful_path{
        my $in=shift;
        if($in !~ /^\//){
                my $t=`pwd`;
                chomp($t);
                return "$t/$in";
        }else{
                return $in;
        }
}
