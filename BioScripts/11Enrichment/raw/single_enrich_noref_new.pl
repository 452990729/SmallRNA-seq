#Edit by jiangxiaoxue
##2012/11/06
## wangshaobin 2013/11/26
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
	-t <s>	miRNA-target-pairs file
	-length <s>	gene length file, which the first col is gene names, and the second is length data
	-goann <s>	gene to GO annotations file
	-koann <s>	gene to KO annotations file
	-o <s>  outdir, eg: sample name	
	-h|?|help	show this help
	-pa <s>         plant or animal
======================================================================================================================================
USAGE
}
my ($help,$input,$target,$length,$goann,$koann,$head,$out,$pa);
GetOptions("h|?|help" => \$help,
	"t=s"=>\$target,
	"length=s"=>\$length,
	"goann=s"=>\$goann,
	"koann=s"=>\$koann,
	"o=s"=>\$out,
	"pa=s"=>\$pa,
);

if(defined($help)||!defined($target)||!defined($length)||!($goann)||!defined($koann)){
	&usage;
	exit(0);
}


my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $pythonExec = $Config->{srnaenv}->{python_v276};

my $add_id="$Bin/add_id.pl";
my $identify="$Bin/identify.py";
my $extract="$Bin/extract_top20.pl";
my $pdf_plot="$Bin/Pathwayscatter_pdf.R";
my $png_plot="$Bin/Pathwayscatter_png.R";

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
open(OUT,">$pwd/$out\_run_GO_enrich.sh");
open(OUT_KEGG,">$pwd/$out\_run_KEGG_enrich.sh");

my $goann_path=get_ful_path($goann);
my $koann_path=get_ful_path($koann);
my $length_path=get_ful_path($length);
$length=basename($length);
$goann=basename($goann);
$koann=basename($koann);

`mkdir -p $pwd/$out/`;
`awk '{print \$2}' $target|sort -u >$pwd/$out/$out.miRNA-geneid`;
`ln -sf $goann_path $pwd/$out/$goann`;
`ln -sf $koann_path $pwd/$out/$koann`;
`ln -sf $length_path $pwd/$out/$length`;

print OUT "echo Start Time:\ndate\n";
print OUT "cd $pwd/$out/\n";

print OUT "echo ================== Run GO enrichment==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "$perlExec $Bin/goseq_graph_v3.pl -i $pwd/$out/$out.miRNA-geneid -goann $pwd/$out/$goann -length $pwd/$out/$length -o $pwd/$out/GO_enrichment/ -p $out\n";

print OUT "<<\\Mark\n";
print OUT "echo ================== Run GO classification==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "sort -k 1,1 $pwd/$out/$goann|join $pwd/$out/$out.miRNA-geneid - |sed 's/ /\\t/g' >$pwd/$out/$out.miRNA-geneid.go.txt\n";
print OUT "mkdir -p GO_classify\n";
print OUT "$perlExec $Bin/go_classification_v3.pl $pwd/$out/$out.miRNA-geneid.go.txt $pwd/$out/GO_classify\n";
print OUT "Mark\n";

print OUT_KEGG "echo ================== Run KEGG enrichment==================\n";
print OUT_KEGG "echo Start Time:\ndate\n";
print OUT_KEGG "mkdir -p $pwd/$out/Pathway/\n";
print OUT_KEGG "cd $pwd/$out/Pathway/\n";
print OUT_KEGG "$perlExec $Bin/runKEGG_enrich_v1.1.pl -diff $pwd/$out/$out.miRNA-geneid -ko $pwd/$out/$koann -g $out -t $pa\n";
print OUT_KEGG "$pythonExec $Bin/pathway_annotation_flow_parallel_simple_tolerant.py --table ${out}.DEG_KEGG_pathway_enrichment_add.xls --abbr ko\n ";
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
