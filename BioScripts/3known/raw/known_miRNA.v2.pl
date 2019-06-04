#Edit by jiangxiaoxue
#2012/09/06
#version2,2013.07.09,by zhangyu

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
example:
perl $0 -r all.map.collapse.fa -f mir_list.txt -o ssc -g 0 -T 10
qsub -cwd -V -l vf=1g -l p=10 ssc_runknown_miRNA.sh

options:
[mandatory parameters]
	-h|?|--help			help information
	-r reads.fa			Input fasta in miRDeep2 format
					eg:	>PAN_123456_x969696
						ATACAATCTACTGTCTTTCCT
	-f <s>				mir_list.txt,a one column list of abbr for species from /PUBLIC/database/RNA/miRBase/miRBase21/organisms.txt which is the same as pathway
[optional parameters]
	-d					permit the mirdeep2 pdfs output, add -d will forbid creat pdf
	-o prj				project name
	-g [int]			number of allowed mismatches when mapping reads to precursors, default 0
	-T [int]			number of alignment threads to launch during bowtie mapping(default:10)
======================================================================================================================================
USAGE
}

my ($help,$d,$reads,$txt,$prj,$mismatch,$threads);
GetOptions(
        "h|?|help"=>\$help,
	"d"=>\$d,
	"r=s"=>\$reads,
	"f=s"=>\$txt,
	"o=s"=>\$prj,
	"g=i"=>\$mismatch,
	"T=i"=>\$threads,
);


my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $mirdeep2 = $Config->{software}->{mirdeep2};
my $srnatoolscli = $Config->{software}->{srnatoolscli};
my $miREvo = $Config->{software}->{miREvo};
my $miRBase21 = $Config->{database}->{miRBase21};
my $parse="$Bin/parse_prefix2fa.pl";
my $ViennaRNA2 = $Config->{software}->{ViennaRNA2};

$prj ||= "prj";
$mismatch = 0 if(!defined($mismatch));
$threads ||= 10;

if(!defined($reads) || !defined($txt) || defined($help)){
	&usage;
	exit 0;
}

if(!defined($d)){
	$d="-d";
}else{
	$d="";
}

my $pwd = `pwd`;
$pwd =~ s/\s+//g;

my $hp="ref_hairpin.fa";
my $mat="ref_mature.fa";

open OUT, ">$prj\_runknown_miRNA.sh" or die $!;
print OUT "echo Start Time:\ndate\n";
print OUT "cd $pwd\n";
#=======================================================edit by jc in 2015/5/16
open TXT,$txt;
while (<TXT>){
	chomp;
	if ($_ eq "all"){
		print OUT "ln -sf $miRBase21/hairpin.fa $hp\n";
		print OUT "ln -sf $miRBase21/mature.fa $mat\n";
	}elsif( $_ eq "animal"){
		print OUT "ln -sf $miRBase21/hairpin_animal.fa $hp\n";
		print OUT "ln -sf $miRBase21/mature_animal.fa $mat\n";
	}elsif( $_ eq "plant"){
		print OUT "ln -sf $miRBase21/hairpin_plant.fa $hp\n";
		print OUT "ln -sf $miRBase21/mature_plant.fa $mat\n";
	}else{
		print OUT "$perlExec $parse $txt $miRBase21/hairpin.fa >$hp\n";
		print OUT "$perlExec $parse $txt $miRBase21/mature.fa >$mat\n";
	}
}
close TXT;
#===========================================================edit by jc in 2015/5/16
if(!(-e "$pwd/$prj")){
	print OUT "mkdir -p $pwd/$prj\n";
}

if(!(-e "$pwd/$prj/ref")){
	print OUT "mkdir -p $pwd/$prj/ref\n";
}
print OUT "echo ================== Format precursor and mature sequence ==================\n";
print OUT "$perlExec $Bin/known_hp_mat_remove_dup.pl $hp $mat $pwd/$prj/ref/known_hairpin_ref.fa $pwd/$prj/ref/known_mature_ref.fa\n\n";

print OUT "echo ================== Run quantifier_gb scan known miRNA ==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "export PATH=$mirdeep2:\$PATH\n";
my $reads_name = basename($reads);
print OUT "ln -sf ",get_ful_path($reads)," $pwd/$prj/$reads_name\n";
print OUT "cd $pwd/$prj/\n";
print OUT "$perlExec $Bin/quantifier_gb_v2.pl -p ref/known_hairpin_ref.fa -m ref/known_mature_ref.fa  -r $reads_name -y $prj.known -g $mismatch -T $threads \n\n";

print OUT "echo ================== Run csv2pdf known miRNA ==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "export PATH=$srnatoolscli:$ViennaRNA2/bin:\$PATH\n";
print OUT "$perlExec $Bin/csv2pdf_gb.pl miRNAs_expressed_all_samples_$prj.known.csv ref/known_hairpin_ref.fa ref/known_mature_ref.fa $prj.known\n";
print OUT "rm $prj.known/image/legend.txt\n";
print OUT "awk '{if(\$4>0){print \$1\"\\t\"\$3}}' miRNAs_expressed_all_samples_$prj.known.csv >$prj.known/hairpin_mature.pairs\n";
print OUT "head -1 $prj.known/miRNAs_expressed_all_samples_$prj.known.csv|awk \'{num=(NF-4)/2;printf(\"miRNA\");for(i=1;i<=num;i++){printf(\"\\t\"\$(i+4))};printf(\"\\n\");}\' >$prj.known/mature.readcount\n";
print OUT "awk \'{if(NR>1){num=(NF-4)/2;printf(\$1\"\\t\"\$2);for(i=1;i<=num;i++){printf(\"\\t\"\$(i+4))};printf(\"\\n\")}}\' $prj.known/miRNAs_expressed_all_samples_$prj.known.csv|sort -k 1,1 -k 2nr,2|awk \'{if(\$1!=name){num=NF-2;printf(\$1);for(i=1;i<=num;i++){printf(\"\\t\"\$(i+2));}printf(\"\\n\");name=\$1}}\' >>$prj.known/mature.readcount\n";
print OUT "awk \'{if(/^>/){title=\$1}else{if(title){if(/total read count/ && \$(NF)>0){print title;title=\"\";marker=1}else{marker=0}}if(marker){print}}}' expression_analyses/expression_analyses_$prj.known/miRBase.mrd >$prj.known/miRBase.mrd\n\n";

print OUT "echo ================== Run Stat known miRNA ==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "$perlExec $Bin/genebwt12count.pl -i expression_analyses/expression_analyses_$prj.known/*_mapped.bwt.ka -r expression_analyses/expression_analyses_$prj.known/precursor.converted -t known_miRNA -o $prj.known -u -s -W\n";
print OUT "awk \'{if(NR==1){print \"Types\\t\"\$0}else if(NR==2){print \"Mapped mature\\t\"\$0}}\' $prj.known/known_miRNA.mapmat.stat >$prj.known/known_miRNA.map.stat\n";

#edit by JC on 2015/12/23================修复known_miRNA.mapref.stat检测到的total前体与 hairpin.fa中学列相等，并且每个样本检测到的前体小于total
print OUT "$perlExec $Bin/stat_known_miRNA_pre.pl $prj.known/miRBase.mrd >$prj.known/known_miRNA.mapref.stat\n";
#=======================================

print OUT "awk \'{if(NR==2){print \"Mapped hairpin\\t\"\$0}}\' $prj.known/known_miRNA.mapref.stat >>$prj.known/known_miRNA.map.stat\n";

print OUT "awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=int(\$i+0.5)}printf(\"Mapped uniq sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"int(\$i+0.5))}printf(\"\\n\");}}\' $prj.known/known_miRNA.uc.stat >>$prj.known/known_miRNA.map.stat\n";
print OUT "awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=int(\$i+0.5)}printf(\"Mapped total sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"int(\$i+0.5))}printf(\"\\n\");}}\' $prj.known/known_miRNA.rc.stat >>$prj.known/known_miRNA.map.stat\n";

## 2014-2-10 add  for unifying "hairpin.fa" and last "*map.stat"
print OUT "cd $prj.known\n";
print OUT "cd ../\n";

print OUT "cp expression_analyses/expression_analyses_$prj.known/*.unmap.fas $prj.known/known_miRNA.unmap.fas\n";
print OUT "$perlExec $Bin/bwt12collapse.pl expression_analyses/expression_analyses_$prj.known/*_mapped.bwt.k1 >$prj.known/known_miRNA.map.fas\n";
print OUT "$perlExec $Bin/bwt_base_bias_plot.pl -i expression_analyses/expression_analyses_$prj.known/*_mapped.bwt.k1\n";
print OUT "cd $pwd\n";
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
