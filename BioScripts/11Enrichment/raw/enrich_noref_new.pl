#Edit by jiangxiaoxue
##2012/11/06
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
example: perl $0 -i DE_dir -c group1vsgroup2,group2vsgroup3 -t miRNA-target.pairs -length geneINFO -goann prj.go.txt -koann koID.annotation -o prj
nohup sh prj_enrich.sh >prj_enrich.sh.o 2>prj_enrich.sh.e &
======================================================================================================================================
Usage: perl $0 [options]
	-i <s>	diff analysis result directoty
	-c <s>	compare directory name, split by ","
	-t <s>	miRNA-target-pairs file
	-length <s>	gene length file, which the first col is gene names, and the second is length data
	-goann <s>	gene to GO annotations file
	-koann <s>	gene to KO annotations file
	-o <s>		project name
	-pa <s>		plant or animal
	-h|?|help	show this help
======================================================================================================================================
USAGE
}


my ($help,$indir,$compares,$target,$length,$goann,$koann,$prj,$t);
GetOptions(
	"h|?|help" => \$help,
	"i=s"=>\$indir,
	"c=s"=>\$compares,
	"t=s"=>\$target,
	"length=s"=>\$length,
	"goann=s"=>\$goann,
	"koann=s"=>\$koann,
	"o=s"=>\$prj,
	"pa=s"=>\$t,
);

if(defined($help)||!defined($indir)||!defined($compares)||!defined($target)||!defined($length)||!($goann)||!defined($koann)||!defined($prj)){
	&usage;
	exit(0);
}


my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};

my $add_id="$Bin/add_id.pl";
my $identify="$Bin/identify.py";
my $extract="$Bin/extract_top20.pl";
my $pdf_plot="$Bin/Pathwayscatter_pdf.R";
my $png_plot="$Bin/Pathwayscatter_png.R";

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
$target=get_ful_path($target);
$indir=get_ful_path($indir);
$length=get_ful_path($length);
$goann=get_ful_path($goann);
$koann=get_ful_path($koann);

open(OUT,">$pwd/$prj\_enrich.sh");
if(!(-e "$pwd/$prj")){
	print OUT "mkdir -p $pwd/$prj/\n";
}
print OUT "cd $pwd/$prj/\n";
print OUT "echo Start Time:\ndate\n";

my @compare=split(",",$compares);
foreach my $i(@compare){
	print OUT "$perlExec $Bin/run_enrich_noref_new.pl -i $indir/$i/$i.DElist.txt -t $target -length $length -goann $goann -koann $koann -o $i -pa $t\n";
	print OUT "qsub -V -cwd -l vf=2G $i\_GO_runenrich.sh\n";
	print OUT "qsub -V -cwd -l vf=2G $i\_KEGG_runenrich.sh\n";
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
