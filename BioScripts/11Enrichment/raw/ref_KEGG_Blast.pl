#Edit by jiangxiaoxue
##2012/11/30
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Config::Tiny;
#my $Bin="/PUBLIC/source/RNA/smallRNA/version3/11Enrichment";
sub usage{
	print STDERR <<USAGE;
======================================================================================================================================
example: perl $0 -i DE_dir -c group1vsgroup2,group2vsgroup3 -t miRNA-target.pairs -species osa -genome genome.fa -gtf genome.gtf -goann gene.go -o prj
nohup sh prj_enrich.sh >prj_enrich.sh.o 2>prj_enrich.sh.e &
======================================================================================================================================
Usage: perl $0 [options]
	-i <s>  diff analysis result directoty
	-t <s>	miRNA-target-pairs file
	-species <s>	abbr for species(from /PUBLIC/software/public/Annotation/kobas2.0-20120208/scripts/sp.list)
	-genome	<s>	genome sequece
	-gtf <s>	Gene sets for your organism, must be a gtf file.
	-genelength <s> Gene length file
	-geneseq <s>      dirGeneSeq
	-o <s>		project name
	-h|?|help	show this help
======================================================================================================================================
USAGE
}
my ($help,$indir,$target,$gtf,$genelength,$dirGeneSeq,$species,$genome,$prj);
GetOptions(
	"h|?|help" => \$help,
	"i=s"=>\$indir,
	"t=s"=>\$target,
	"gtf=s"=>\$gtf,
	"genelength=s"=>\$genelength,
	"geneseq=s"=>\$dirGeneSeq,
	"species=s"=>\$species,
	"genome=s"=>\$genome,
	"o=s"=>\$prj,
);

if(defined($help)||!defined($indir)||!defined($target)||!defined($species)||!defined($prj)){
	&usage;
	exit(0);
}


my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $kobas_v20140801 = $Config->{database}->{kobas_v20140801};
my $kobas = $Config->{database}->{kobas};

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
$target=get_ful_path($target);
$indir=get_ful_path($indir);
if(defined $genome){$genome=get_ful_path($genome);}
if(defined $gtf){$gtf=get_ful_path($gtf);}
if(defined $genelength){$genelength=get_ful_path($genelength);}

open(OUT,">$pwd/run_Balst.sh");

print OUT "echo Start Time:\ndate\n";

if(defined $gtf){
                print OUT"PYTHONPATH=$kobas_v20140801/src:$kobas:\$PYTHONPATH\n";
		print OUT "$perlExec $Bin/run_ref_KEGG_Blast.pl -i $indir -t $target -species $species -genome $genome -gtf $gtf -o KEGG\n";
		print OUT "sh KEGG_Blast.sh\n";
	}
if(defined $genelength){
		print OUT "$perlExec $Bin/run_enrich_ref_new.pl -i $indir -t $target -species $species  -length $genelength -geneseq $dirGeneSeq -o KEGG\n";
                print OUT "sh KEGG_Blast.sh\n";
	}
print OUT "echo End Time:\ndate\n";
print OUT "cd $pwd\n";
close(OUT);

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
