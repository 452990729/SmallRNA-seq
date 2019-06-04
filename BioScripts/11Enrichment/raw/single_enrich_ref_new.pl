#edit by zhangyu,2013.07.05
#wangshaobin 2013.11.8
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
        -t <s>  miRNA-target-pairs file
        --head  gene-list with headinfo
        -species <s>    abbr for species(from /PUBLIC/software/public/Annotation/kobas2.0-20120208/scripts/sp.list)
        -genome <s>     genome sequece
        -gtf <s>        Gene sets for your organism, must be a gtf file.
        -goann <s>      gene to GO annotations file
        -o <s>  project output directory
        -h|?|help       show this help
======================================================================================================================================
USAGE
}
my ($help,$target,$gtf,$species,$goann,$genome,$head,$pro);
GetOptions("h|?|help" => \$help,
        "t=s"=>\$target,
        "gtf=s"=>\$gtf,
        "species=s"=>\$species,
        "genome=s"=>\$genome,
        "goann=s"=>\$goann,
        "head"=>\$head,
        "o=s"=>\$pro,
);

if(defined($help)||!defined($target)||!defined($gtf)||!defined($species)||!defined($goann)||!defined($genome)||!defined($pro)){
        &usage;
        exit(0);
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $pythonExec = $Config->{srnaenv}->{python_v276};
my $pwd = `pwd`;
$pwd =~ s/\s+//g;

open(OUT,">$pwd/$species\_run_GO_enrich.sh");
open(OUT_KEGG,">$pwd/$species\_run_KEGG_enrich.sh");

my $goann_path=get_ful_path($goann);
my $genome_path=get_ful_path($genome);
my $gtf_path=get_ful_path($gtf);

$gtf=basename($gtf);
$goann=basename($goann);
$genome=basename($genome);

`mkdir -p $pwd/$pro/`;
#print OUT "mkdir $pwd/$pro/\n";
print OUT "echo Start Time:\ndate\n";
#print OUT "echo ================== preparing data==================\n";
`awk '{print \$2}' $target | sort -u >$pwd/$pro/$pro.miRNA-geneid`;
`ln -sf $goann_path $pwd/$pro/$goann`;
`ln -sf $genome_path $pwd/$pro/$genome`;
`ln -sf $gtf_path $pwd/$pro/$gtf`;

print OUT "cd $pwd/$pro/\n";

print OUT "echo ================== Run GO enrichment==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "$perlExec $Bin/goseq_graph.pl -i $pwd/$pro/$pro.miRNA-geneid -goann $pwd/$pro/$goann -gtf $pwd/$pro/$gtf -o $pwd/$pro/GO_enrichment -p $pro\n";
print OUT "echo ================== Run GO classification==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "sort -k 1,1 $pwd/$pro/$goann|join $pwd/$pro/$pro.miRNA-geneid - |sed 's/ /\\t/g' >$pwd/$pro/$pro.miRNA-geneid.go.txt\n";
print OUT "mkdir -p $pwd/$pro/GO_classify\n";
print OUT "$perlExec $Bin/go_classification_v3.pl $pwd/$pro/$pro.miRNA-geneid.go.txt $pwd/$pro/GO_classify\n";

print OUT_KEGG "echo ================== Run KEGG===========================\n";
print OUT_KEGG "echo Start Time:\ndate\n";
print OUT_KEGG "mkdir -p $pwd/$pro/Pathway/\n";
print OUT_KEGG "cd $pwd/$pro/\n";
print OUT_KEGG "$perlExec $Bin/runKEGG.pl -directory-gtf $gtf -gene-id $pro.miRNA-geneid -exist-genome y -directory-kobas-output $pwd/$pro/Pathway/ -directory-genome $genome -species $species -samplename $pro\n";
print OUT_KEGG "cd $pwd/$pro/Pathway/\n";
print OUT_KEGG "$pythonExec $Bin/pathway_annotation_flow_parallel_simple_tolerant.py --table add.${pro}.identify.xls --abbr $species\n ";

print OUT_KEGG "echo End Time:\ndate\n";

close(OUT);
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

