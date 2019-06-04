#Edit by jiangxiaoxue
##2012/11/06
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
	-length <s>	gene length file, which the first col is gene names, and the second is length data
	-goann <s>	gene to GO annotations file
	-koann <s>	gene to KO annotations file
	-o <s>	project name(suggest using the different express analysis each group name like sC1_vs_sC2)
	-pa <s> 		plant or animal
	-h|?|help	show this help
======================================================================================================================================
USAGE
}
my ($help,$input,$target,$length,$goann,$koann,$head,$com,$t);
GetOptions("h|?|help" => \$help,
	"i=s"=>\$input,
	"t=s"=>\$target,
	"length=s"=>\$length,
	"goann=s"=>\$goann,
	"koann=s"=>\$koann,
	"head"=>\$head,
	"o=s"=>\$com,
	"pa=s"=>\$t,
);

if(defined($help)||!defined($input)||!defined($target)||!defined($length)||!($goann)||!defined($koann)){
	&usage;
	exit(0);
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $R_v2153 = $Config->{srnaenv}->{R_v2153};
my $LIBRARY_PATH = $Config->{srnaenv}->{LIBRARY_PATH};
my $pythonExec = $Config->{srnaenv}->{python_v2710};
my $Site_Packages = $Config->{srnaenv}->{Site_Packages};
my $kobas_v20120208 = $Config->{database}->{kobas_v20120208};
my $kobas_v20140801 = $Config->{database}->{kobas_v20140801};
my $kobas = $Config->{database}->{kobas};

my $add_id="$Bin/add_id.pl";
my $identify="$Bin/identify.py";
my $extract="$Bin/extract_top20.pl";
my $pdf_plot="$Bin/Pathwayscatter_pdf.R";
my $png_plot="$Bin/Pathwayscatter_png.R";
my $goann_path=get_ful_path($goann);
my $koann_path=get_ful_path($koann);
my $length_path=get_ful_path($length);
$length=basename($length);
$goann=basename($goann);
$koann=basename($koann);


my $pwd = `pwd`;
$pwd =~ s/\s+//g;
open(OUT,">$pwd/$com\_GO_runenrich.sh");
open(OUT_KEGG,">$pwd/$com\_KEGG_runenrich.sh");

`mkdir -p $pwd/$com/`;
if(defined($head)){
	`awk '{if(NR>1){print}}' $input |sort -k 1,1 >$pwd/$com/$com.diffmiRNAID`;
}else{
	`sort -k 1,1 $input >$pwd/$com/$com.diffmiRNAID`;
}
`sort -k 1,1 $target|join $pwd/$com/$com.diffmiRNAID - |sort -u >$pwd/$com/$com.diffmiRNA-gene.pairs`;
`awk '{print \$2}' $pwd/$com/$com.diffmiRNA-gene.pairs|sort -u >$pwd/$com/$com.diffmiRNA-geneid`;
`ln -sf $goann_path $pwd/$com/$goann`;
`ln -sf $koann_path $pwd/$com/$koann`;
`ln -sf $length_path $pwd/$com/$length`;

print OUT "echo ================== Run different express miRNA target gene enrichment analysis==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "cd $pwd/$com/\n";

print OUT "echo ================== Run GO enrichment==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "$perlExec $Bin/goseq_graph_v3.pl -i $pwd/$com/$com.diffmiRNA-geneid -goann $pwd/$com/$goann -length $pwd/$com/$length -o $pwd/$com/GO2/ -p $com\n";


print OUT "echo ================== Run GO classification==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "sort -k 1,1 $pwd/$com/$goann|join $pwd/$com/$com.diffmiRNA-geneid - |sed 's/ /\\t/g' >$pwd/$com/$com.diffmiRNA-geneid.go.txt\n";
print OUT "mkdir -p GO_classify\n";
print OUT "$perlExec $Bin/go_classification_v3.pl $pwd/$com/$com.diffmiRNA-geneid.go.txt $pwd/$com/GO_classify\n";

print OUT_KEGG "echo ================== Run KEGG enrichment==================\n";
print OUT_KEGG "echo Start Time:\ndate\n";
print OUT_KEGG "mkdir -p $pwd/$com/Pathway/\n";
print OUT_KEGG "export PATH=$R_v2153:\$PATH\n";
print OUT_KEGG "export LD_LIBRARY_PATH=$LIBRARY_PATH/R_v2153:\$LD_LIBRARY_PATH\n";
print OUT_KEGG "export PYTHONPATH=$Site_Packages:$kobas_v20140801/src:$kobas_v20120208/src:$kobas:\$PYTHONPATH\n";
print OUT_KEGG "cd $pwd/$com/Pathway/\n";
print OUT_KEGG "$perlExec $Bin/runKEGG_enrich_v1.1.pl -diff $pwd/$com/$com.diffmiRNA-geneid -ko $pwd/$com/$koann -g $com -t $t\n";
print OUT_KEGG "$pythonExec $Bin/pathway_annotation_flow_parallel_simple_tolerant.py  --table $pwd/$com/Pathway/$com.DEG_KEGG_pathway_enrichment_add.xls --abbr ko\n";#edit in 2015/3/12
=head
print OUT "echo ================== Run KEGG enrichment==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "mkdir $pwd/$com/Pathway/\n";
print OUT "perl $Bin/join.diffgeneid.annot.pl $pwd/$com/$koann $pwd/$com/$com.diffmiRNA-geneid $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation\n";
print OUT "python $identify -f $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation -b $pwd/$com/$koann -o $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify\n";
#print OUT "python $identify -f $pwd/$com/$koann $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation -d K -m h -n BH -b ko -s not standard -o $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify\n";
print OUT "perl $add_id $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify >$pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add\n";
print OUT "perl $extract $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add >$pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add.top20\n";
print OUT "Rscript $pdf_plot $com $pwd/$com/Pathway/ $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add.top20\n";
print OUT "Rscript $png_plot $com $pwd/$com/Pathway/ $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add.top20\n";
print OUT "cp $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add.top20 $pwd/$com/Pathway/$com.Pathscatter.txt\n";
print OUT "cd $pwd/$com/Pathway/\n";
print OUT "python /PUBLIC/source/RNA/RefRNA/DGE/Enrichment/pathway_annotation_flow_parallel_annotationfault_tolerant.pyc --table $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add\n";
#print OUT "cut -f 1-9,12-13 $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add >$pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add.cut\n";
#print OUT "python $Bin/pathway_annotation_flow.py $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add.cut $Bin\n";
#print OUT "cp $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add.cut $pwd/$com/Pathway/$com.pathway.xls\n";
#print OUT "mv $pwd/$com/$com.diffmiRNA-gene.koID.annotation.identify.add.cut_rendered_html_new.html $pwd/$com/Pathway/$com.pathway.html\n";
#print OUT "cp -r $Bin/src $pwd/$com/Pathway/\n";
#print OUT "python $Bin/pathway_annotation_flow.py $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add\n";
#print OUT "cp $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add $pwd/$com/Pathway/$com.pathway.xls\n";
#print OUT "mv $pwd/$com/Pathway/$com.diffmiRNA-gene.koID.annotation.identify.add_rendered_html_new.html $pwd/$com/Pathway/$com.pathway.html\n";
=cut
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
