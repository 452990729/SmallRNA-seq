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
	-goann <s>	gene to GO annotations file
	-ann <s>        add description
	-o <s>	project name(suggest using the different express analysis each group name like sC1_vs_sC2)
	-h|?|help	show this help
======================================================================================================================================
USAGE
}
my ($help,$input,$target,$gtf,$species,$goann,$genome,$head,$com,$ann);
GetOptions("h|?|help" => \$help,
	"i=s"=>\$input,
	"t=s"=>\$target,
	"gtf=s"=>\$gtf,
	"species=s"=>\$species,
	"genome=s"=>\$genome,
	"goann=s"=>\$goann,
	"head"=>\$head,
	"ann=s"=>\$ann,
	"o=s"=>\$com,
);

if(defined($help)||!defined($input)||!defined($target)||!defined($gtf)||!defined($species)||!defined($goann)||!defined($genome)){
	&usage;
	exit(0);
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $pythonExec = $Config->{srnaenv}->{python_v2710};
my $R_v2153 = $Config->{srnaenv}->{R_v2153};
my $LIBRARY_PATH = $Config->{srnaenv}->{LIBRARY_PATH};
my $Site_Packages = $Config->{srnaenv}->{Site_Packages};
my $kobas_v20120208 = $Config->{database}->{kobas_v20120208};
my $kobas_v20140801 = $Config->{database}->{kobas_v20140801};
my $kobas = $Config->{database}->{kobas};


my $goann_path=get_ful_path($goann);
my $genome_path=get_ful_path($genome);
my $gtf_path=get_ful_path($gtf);

$goann=basename($goann);
$genome=basename($genome);
$gtf=basename($gtf);

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
`$perlExec $Bin/add_description.pl $pwd/$com/$com.diffmiRNA-gene.pairs $ann`;
`mv $pwd/$com/$com.diffmiRNA-gene.pairs.mg $pwd/$com/$com.diffmiRNA-gene.pairs`;
`awk '{if(NR>0){print \$2}}' $pwd/$com/$com.diffmiRNA-gene.pairs|sort -u >$pwd/$com/$com.diffmiRNA-geneid`;
`ln -sf $goann_path $pwd/$com/$goann`;
`ln -sf $genome_path $pwd/$com/$genome`;
`ln -sf $gtf_path $pwd/$com/$gtf`;
print OUT "echo ================== Run different express miRNA target gene enrichment analysis==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "cd $pwd/$com/\n";

print OUT "echo ================== Run GO enrichment==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "$perlExec $Bin/goseq_graph_v3.pl  -i $pwd/$com/$com.diffmiRNA-geneid -goann $pwd/$com/$goann -gtf $pwd/$com/$gtf -o $pwd/$com/GO2/ -p $com\n";

=head edit by jc on 2016/1/5
print OUT "echo ================== Run GO classification==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "sort -k 1,1 $pwd/$com/$goann|join $pwd/$com/$com.diffmiRNA-geneid - |sed 's/ /\\t/g' >$pwd/$com/$com.diffmiRNA-geneid.go.txt\n";
print OUT "mkdir $pwd/$com/GO_classify\n";
print OUT "perl /PUBLIC/software/RNA/gene_annotation/scripts/go_classification_v3.pl $pwd/$com/$com.diffmiRNA-geneid.go.txt $pwd/$com/GO_classify\n";
=cut
print OUT_KEGG "echo ================== Run KOBAS==================\n";
print OUT_KEGG "echo Start Time:\ndate\n";
print OUT_KEGG "mkdir -p $pwd/$com/Pathway/\n";
print OUT_KEGG "export PATH=$R_v2153:\$PATH\n";
print OUT_KEGG "export LD_LIBRARY_PATH=$LIBRARY_PATH/R_v2153:\$LD_LIBRARY_PATH\n";
print OUT_KEGG "export PYTHONPATH=$Site_Packages:$kobas_v20140801/src:$kobas_v20120208/src:$kobas:\$PYTHONPATH\n";
print OUT_KEGG "cd $pwd/$com/\n";
print OUT_KEGG "$perlExec $Bin/runKEGG_v1.pl  -directory-gtf $gtf -gene-id $com.diffmiRNA-geneid -exist-genome y -directory-kobas-output $pwd/$com/Pathway/ -directory-genome $genome -species $species -samplename $com\n";
print OUT_KEGG "cd $pwd/$com/Pathway/\n";
print OUT_KEGG "$pythonExec $Bin/pathway_annotation_flow_parallel_simple_tolerant.py --table add.$com.identify.xls --abbr $species\n";
print OUT_KEGG "echo End Time:\ndate\n";
print OUT_KEGG "cd $pwd\n";
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
