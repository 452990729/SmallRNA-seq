#Edit by jiangxiaoxue
##2012/10/12
#modify by wangshaobin

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Config::Tiny;


sub usage{
        print STDERR <<USAGE;
======================================================================================================================================
example:
perl $0 -m known_mature.fa,novel_mature.fa -r 3UTR.fa -o ssr
qsub -cwd -V -l vf=1g ssr_runtarget.sh
=========================================================================
Usage: perl $0 [options]
options:
	-h|?|--help	help information
	-m	known_mature.fa,novel_mature.fa(single or mutiple split by "," is ok)
	-r	3UTR.fa for animals, the three_prime_UTR fasta
	-a	a two raw gene annotation file, e.g,geneID\tdescription
	-gtf    gtf
	-o	prj
=========================================================================
USAGE
}


my ($help,$type,$refer,$gtf,$matures,$anno,$prj);
GetOptions(
	"h|?|help"=>\$help,
	"r=s"=>\$refer,
	"gtf=s"=>\$gtf,
	"m=s"=>\$matures,
	"type=s"=>\$type,
	"a=s"=>\$anno,
	"o=s"=>\$prj,
);

if(!defined($refer) || !defined($gtf) || !defined($type) || !defined($matures) || !defined($prj) ||!defined($anno) ||defined($help)){
	&usage;
	exit 0;
}


my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $pythonExec = $Config->{srnaenv}->{python_v276};
my $LIBRARY_PATH = $Config->{srnaenv}->{LIBRARY_PATH};
my $miranda_path="$Bin/runtarget_split1.py";
#my $RNAhybrid_path="/PROJ/RNA/share/software/RNAhybrid/bin/RNAhybrid";


my $pwd = `pwd`;
my $dir =`dir`;
$pwd =~ s/\s+//g;

$matures=~s/,/ /g;

#my  $indir=get_ful_path($refer);

open(OUT,">$pwd/$prj\_runtarget.sh");
print OUT "#!/bin/sh \n";
print OUT "export LD_LIBRARY_PATH=$LIBRARY_PATH/HTSeq_Atlas:\$LD_LIBRARY_PATH \n";
print OUT "dir=$pwd/\n";
print OUT "mkdir -p $pwd/$prj\n";
print OUT "ln -sf ",get_ful_path($refer)," $prj/\n";
$refer=basename($refer);
print OUT "cat ",$matures," >$prj/$prj\_mature.fa\n";
print OUT "cd $pwd/$prj\n";
print OUT "$pythonExec $miranda_path --mature $prj\_mature.fa --refRNA $refer --type $type --outdir $pwd/$prj \n";
print OUT "$perlExec $Bin/trans_gene_pairs.pl -gtf $gtf -tp Common/all_target.xls -o Common/all_target_gene.xls\n";
print OUT "head -15 Common/commom_target.xls >Common/commom_target_example.xls\n";
print OUT "awk '{print \$1\"\\t\"\$3}' Common/all_target_gene.xls|sort -u >Common/all_targets_gene.pairs\n";
print OUT "cp Common/all_targets_gene.pairs ./mature.miranda_targets_gene.pairs \n";

print OUT "$perlExec $Bin/targets.pairs_transfer_v4.pl Common/all_targets_gene.pairs  $anno Common/all_target_gene.xls\n";
print OUT "$pythonExec $Bin/target_clean.py --outdir $pwd/$prj --abbr $prj  --version v2.3 \n";

close(OUT);

sub get_ful_path{
        my $in=shift;
        if($in !~ /^\//){
                my $t=`pwd`;
                chomp($t);
                return "$t/$in";
        }
        else{
                return $in;
        }
}

