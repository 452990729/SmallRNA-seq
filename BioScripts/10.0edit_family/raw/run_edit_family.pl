#Edit by jiangxiaoxue
#2012/09/25

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Config::Tiny;

sub usage{
	print STDERR <<USAGE;
======================================================================================================================================
example:perl $0 --known_dir ../4k_miRNA/prj/prj.known --novel_dir ../8n_miRNA/prj/prj.novel -i ../3map/prj/SP1.collapse.fa,../3map/prj/SP2.collapse.fa -n SP1,SP2 -o prj
nohup sh prj_runedit_family.sh >prj_runedit_family.sh.o 2>prj_runedit_family.sh.e &

Usage: perl $0 [options]
options:
	-h|?|--help	help information
	--known_dir	the prj.known from the third step known miRNA analysis
	--novel_dir	the prj.novel from the eight step novel miRNA analysis
	-i read1.fa,read2.fa	Input fasta in miRDeep2 format
                                        eg:     >PAN_123456_x969696
                                                ATACAATCTACTGTCTTTCCT
	different file seperated by "," corresponding to the next sample name list
	-n sp1,sp2	sample name list,seperated by ","
	-o outprj	give a project of this analysis samples group
======================================================================================================================================
USAGE
}

my ($help,$known_dir,$novel_dir,$reads,$samples,$prj);
GetOptions(
	"h|?|help"=>\$help,
	"i=s"=>\$reads,
	"n=s"=>\$samples,
	"known_dir=s"=>\$known_dir,
	"novel_dir=s"=>\$novel_dir,
	"o=s"=>\$prj
);

if(defined($help)||!defined($reads)||!defined($samples)||!defined($known_dir)||!defined($novel_dir)||!defined($prj)){
	&usage;
	exit 0;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $python_v276 = $Config->{srnaenv}->{python_v276};
my $python_v2710 = $Config->{srnaenv}->{python_v2710};
my $mirdeep2 = $Config->{software}->{mirdeep2};
my $miRBase21 = $Config->{database}->{miRBase21};
my $Rfam = $Config->{database}->{Rfam};
my $blast = $Config->{software}->{blast1};
my $infernal = $Config->{software}->{infernal};
my $rfam_pl="$Bin/rfam_scan.pl";

my $pwd = `pwd`;
$pwd =~ s/\s+//g;

open OUT, ">$prj\_runedit_family.sh" or die $!;
print OUT "echo Start Time:\ndate\n";
if(!(-e "$pwd/$prj")){
	print OUT "mkdir -p $pwd/$prj\n";
}

print OUT "echo ================== Run Stat $prj edit analysis ==================\n";
print OUT "cat $known_dir/hairpin.fa >$pwd/$prj/$prj\_known\_hairpin.fa\n";
print OUT "cat $novel_dir/hairpin.fa >$pwd/$prj/$prj\_novel\_hairpin.fa\n";
print OUT "cat $known_dir/hairpin.fa $novel_dir/hairpin.fa >$pwd/$prj/$prj\_hairpin.fa\n";
print OUT "cat $known_dir/hairpin_mature.fa $novel_dir/hairpin_mature.fa >$pwd/$prj/$prj\_hairpin_mature.fa\n";
print OUT "cd $pwd/$prj/\n";

print OUT "echo ================== Run Stat $prj edit analysis ==================\n";
print OUT "echo Start Time:\ndate\n";
my @sample=(split /,/, $samples);
my @read=(split /,/, $reads);
print OUT "export PATH=$mirdeep2:\$PATH\n";
print OUT "#export PATH=$mirdeep2/lib:\$PERL5LIB\n";
foreach my $i(0..$#sample){
	print OUT "ln -sf ",get_ful_path($read[$i])," $pwd/$prj/$sample[$i].collapse.fa\n";
	print OUT "$perlExec $Bin/quantifier_gb.pl -p $prj\_hairpin.fa -m $prj\_hairpin_mature.fa -r $sample[$i].collapse.fa -y $sample[$i] -g 2 -T 10 -d\n";
	print OUT "$python_v276 $Bin/miRNA_editing.py  expression_analyses/expression_analyses_$sample[$i]/miRBase.mrd $sample[$i]\n";
	print OUT "awk -F\"\\t\" '{if(/^>/){if(\$2){marker=1}else{marker=0}}if(marker){print}}'  $sample[$i]_editing_stats.txt|head -10 >$sample[$i]_editing_stats.example.txt\n";
}

print OUT "echo ================== Run Stat $prj family analysis ==================\n";
print OUT "echo Start Time:\ndate\n";
#print OUT "PYDIR=\"/PUBLIC/software/public/python-2.7.8\"\n";
#print OUT "export PATH=\${PYDIR}/bin:\$PATH\n";
#print OUT "export PYTHONPATH=\${PYDIR}/lib/python2.7:\$PYTHONPATH\n";
#print OUT "export LD_LIBRARY_PATH=\${PYDIR}/lib:\$LD_LIBRARY_PATH\n";
print OUT "export PATH=$blast:$infernal/src:\$PATH\n";
print OUT "$perlExec $rfam_pl -blastdb $Rfam/Rfam.fasta $Rfam/Rfam.cm $prj\_novel\_hairpin.fa -bt 1e-05 -o novel_hairpin.rfam.gff3\n";
print OUT "$python_v2710 $Bin/miRNA_family.py --known $prj\_known\_hairpin.fa --unknown $prj\_novel\_hairpin.fa --prefix $prj --family_hp $miRBase21/hairpin.fa --family_miFam $miRBase21/miFam.dat --novo_gff novel_hairpin.rfam.gff3 \n";
print OUT "$perlExec $Bin/miRNA_family_show.pl $prj\_miRNA_family.txt $prj\_miRNA_family\n";
print OUT "cut -f 1-10 $prj\_miRNA_family.mir_num.txt|awk -F\"\\t\" '{if(\$2>0 || \$3>0){print}}'|head >$prj\_miRNA_family.mir_num.example.txt\n"; 
print OUT "cut -f 1-10 $prj\_miRNA_family.mir_sign.txt|awk -F\"\\t\" '{if(\$2==\"+\" || \$3==\"+\"){print}else if(NR==1){print}}'|head >$prj\_miRNA_family.mir_sign.example.txt\n";
print OUT "cd $pwd/\n";
print OUT "echo End Time:\ndate\n";


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
