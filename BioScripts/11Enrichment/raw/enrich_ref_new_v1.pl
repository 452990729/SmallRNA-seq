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
example: perl $0 -i DE_dir -c group1vsgroup2,group2vsgroup3 -t miRNA-target.pairs -species osa -genome genome.fa -gtf genome.gtf -goann gene.go -o prj
nohup sh prj_enrich.sh >prj_enrich.sh.o 2>prj_enrich.sh.e &
======================================================================================================================================
Usage: perl $0 [options]
	-i <s>  diff analysis result directoty
	-c <s>  compare directory name, split by ","
	-t <s>	miRNA-target-pairs file
	-species <s>	abbr for species(from /PUBLIC/software/public/Annotation/kobas2.0-20120208/scripts/sp.list)
	-genome	<s>	genome sequece
	-gtf <s>	Gene sets for your organism, must be a gtf file.
	-genelength <s> Gene length file
	-geneseq <s>      dirGeneSeq
	-goann <s>	gene to GO annotations file
	-o <s>		project name
	-ann <s>         genenamefile
	-h|?|help	show this help
======================================================================================================================================
USAGE
}
my ($help,$indir,$compares,$target,$gtf,$genelength,$dirGeneSeq,$species,$goann,$genome,$prj,$ann);
GetOptions(
	"h|?|help" => \$help,
	"i=s"=>\$indir,
	"c=s"=>\$compares,
	"t=s"=>\$target,
	"gtf=s"=>\$gtf,
	"genelength=s"=>\$genelength,
	"geneseq=s"=>\$dirGeneSeq,
	"species=s"=>\$species,
	"genome=s"=>\$genome,
	"goann=s"=>\$goann,
	"ann=s"=>\$ann,
	"o=s"=>\$prj,
);

if(defined($help)||!defined($indir)||!defined($compares)||!defined($target)||!defined($species)||!defined($goann)||!defined($prj)){
	&usage;
	exit(0);
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $pwd = `pwd`;
$pwd =~ s/\s+//g;
$target=get_ful_path($target);
$indir=get_ful_path($indir);
$goann=get_ful_path($goann);
if(defined $genome){$genome=get_ful_path($genome);}
if(defined $gtf){$gtf=get_ful_path($gtf);}
if(defined $genelength){$genelength=get_ful_path($genelength);}
if(defined $ann){$ann=get_ful_path($ann);}


open(OUT,">$pwd/$prj\_enrich.sh");
if(!(-e "$pwd/$prj")){
	 print OUT "mkdir $pwd/$prj/\n";
}
print OUT "cd $pwd/$prj/\n";
print OUT "echo Start Time:\ndate\n";

my @compare=split(",",$compares);
foreach my $i(@compare){
	if(defined $gtf){
		print OUT "$perlExec $Bin/run_enrich_ref_new_v2.pl -i $indir/$i/$i.DElist.txt -t $target -species $species -genome $genome -gtf $gtf -goann $goann -o $i -ann $ann\n";
		print OUT "qsub -V -cwd -l vf=2G $i\_GO_runenrich.sh\n";
		print OUT "qsub -V -cwd -l vf=2G $i\_KEGG_runenrich.sh\n";
	}
	if(defined $genelength){
		print OUT "$perlExec $Bin/run_enrich_ref_new.pl -i $indir/$i/$i.DElist.txt -t $target -species $species  -length $genelength -geneseq $dirGeneSeq -goann $goann -o $i\n";
		print OUT "qsub -V -cwd -l vf=2G $i\_GO_runenrich.sh\n";
                print OUT "qsub -V -cwd -l vf=2G $i\_KEGG_runenrich.sh\n";
	}
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
