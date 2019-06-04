#Edit by jiangxiaoxue
##2012/10/12

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
perl $0 -m known_mature.fa,novel_mature.fa -r Trinity.fasta -f best_candidates.eclipsed_orfs_removed.gff3 -o sCY
qsub -cwd -V -l vf=200M sCY_runtarget.sh
=========================================================================
Usage: perl $0 [options]
options:
	-h|?|--help	help information
	-m	known_mature.fa,novel_mature.fa(single or mutiple split by "," is ok)
	-r	Trinity.fasta
		form animals's transcriptome assembly of trinity
	-f	Trinity_cds_gff3
		form transcriptome Trinity_cds prediction(best_candidates.eclipsed_orfs_removed.gff3)
	-a	Blast_NT.xls
		from transcriptome Annotation
	-o	prj
=========================================================================
USAGE
}

my ($help,$refer,$gff3,$matures,$anno,$prj);
GetOptions(
	"h|?|help"=>\$help,
	"r=s"=>\$refer,
	"f=s"=>\$gff3,
	"m=s"=>\$matures,
	"a=s"=>\$anno,
	"o=s"=>\$prj,
);

if(!defined($refer) || !defined($gff3) || !defined($matures) || !defined($prj) || !defined($anno) ||defined($help)){
	&usage;
	exit 0;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $miranda = $Config->{software}->{miranda};

my $pwd = `pwd`;
$pwd =~ s/\s+//g;

$matures=~s/,/ /g;
=head
my @mature=(split /,/, $matures);
if(@mature != 2){
	die "please input -m as known_mature.fa,novel_mature.fa";
}
=cut

open(OUT,">$pwd/$prj\_runtarget.sh");
print OUT "mkdir -p $pwd/$prj\n";
print OUT "ln -sf ",get_ful_path($refer)," $prj/\n";
print OUT "ln -sf ",get_ful_path($gff3)," $prj/\n";
print OUT "cat ",$matures," >$prj/$prj\_mature.fa\n";
#print OUT "cat ",get_ful_path($mature[0])," ",get_ful_path($mature[1])," >$prj/$prj\_mature.fa\n";
print OUT "cd $pwd/$prj\n";
print OUT "$perlExec $Bin/Trinity_gff3-3UTR.pl ",basename($refer)," ",basename($gff3)," 16 >Trinity_3UTR.fasta\n";
print OUT "$miranda $prj\_mature.fa Trinity_3UTR.fasta -sc 140 -en -10 -scale 4 -strict -out $prj\_mature.miranda_targets_out\n";
print OUT "$perlExec $Bin/mi_result_fmt.pl $prj\_mature.miranda_targets_out $prj\_mature.miranda_targets_out.fmt\n";
print OUT "$perlExec $Bin/mi_result.pl -i $prj\_mature.miranda_targets_out -o $prj\_mature.miranda_targets\n";
#edit by wanghuizhen,2016-5-4 
print OUT "awk 'BEGIN{print \"miRNA\\ttarget_gene\"}{if(name!=\$1){print}}' $prj\_mature.miranda_targets.pairs |head -15 >$prj\_mature.miranda_targets.pairs.example\n";
print OUT "cp $prj\_mature.miranda_targets.pairs.example $prj\_mature.miranda_targets_example.xls\n";
print OUT "cp $prj\_mature.miranda_targets.pairs $prj\_mature.miranda_targets_gene.pairs\n";

print OUT "$perlExec $Bin/targets.pairs_transfer_v2.pl $prj\_mature.miranda_targets.pairs $anno\n";

print OUT "cd $pwd\n";
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

