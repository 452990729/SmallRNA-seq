#Edit by jiangxiaoxue
##2012/09/06
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
for animals:
	perl $0 -i former.unmap.fas -r genome.fas(bowtie index prefix) -m /PUBLIC/database/RNA/miRBase/miRBase20/mature_animal.fa -o ssc
	add -c option can more quick:
		perl $0 -i former.unmap.fas -r genome.fas(bowtie index prefix) -m /PUBLIC/database/RNA/miRBase/miRBase20/mature_animal.fa -o ssc -c
qsub -cwd -V -l vf=200M ssc_runnovel_miRNA.sh
for plants:
	monocot:
		perl $0 -i former.unmap.fas -r genome.fas(bowtie index prefix) -m /PUBLIC/database/RNA/miRBase/miRBase20/mature_plant.fa -o ssc --mode 2	[like animals,add -c option can more quick]
	dicots.
		perl $0 -i former.unmap.fas -r genome.fas(bowtie index prefix) -m /PUBLIC/database/RNA/miRBase/miRBase20/mature_plant.fa -o ssc --mode 3	[like animals,add -c option can more quick]

options:
[mandatory parameters]
	-h|?|--help			help information
	-i reads.fa			Input fasta in miRDeep2 format
					eg:	>PAN_123456_x969696
						ATACAATCTACTGTCTTTCCT
	-r  <str>			The prefix of the bowtie index, for the genome reference file
	-m mature.fa			miRNA sequences from miRBase
[optional parameters]
	--mode <1/2/3>			predicttion model, 1: animal; 2: monocot; 3 dicots.  default=1"
	-k  <int>			a read is allowed to map up to this number of positions in the genome, default=5(if mode==1); for plant, default=15(if mode==2 or 3)
	-c				disable randfold analysis during novel predict(for large genome choose this para, more quick)
	-d				permit the mirdeep2 pdfs output
	-o prj				project name(default prj)
	-v [int]			number of allowed mismatches when mapping reads to precursors, default 0
	-g [int]			maximum number of precursors to analyze when automatic excision gearing is used."
					default=50000, if set to -1 all precursors will be analyzed"
	-T [int]			number of alignment threads to launch during bowtie mapping(default:10)
======================================================================================================================================
USAGE
}

my ($help,$c,$d,$reads,$genome,$mat,$mode,$k,$prj,$mismatch,$gear,$threads);
GetOptions(
        "h|?|help"=>\$help,
	"c"=>\$c,
	"d"=>\$d,
	"i=s"=>\$reads,
	"r=s"=>\$genome,
	"m=s"=>\$mat,
	"mode=i"=>\$mode,
	"k=i"=>\$k,
	"o=s"=>\$prj,
	"v=i"=>\$mismatch,
	"g=i"=>\$gear,
	"T=i"=>\$threads,
);

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $bowtie = $Config->{software}->{bowtie1};
my $mirdeep2 = $Config->{software}->{mirdeep2};
my $srnatoolscli = $Config->{software}->{srnatoolscli};
my $miREvo = $Config->{software}->{miREvo};
my $RNAfold = $Config->{software}->{RNAfold};
my $ViennaRNA2 = $Config->{software}->{ViennaRNA2};


$prj ||= "prj";
$mismatch = 0 if(!defined($mismatch));
$gear= 50000 if(!defined($gear));
$threads ||= 10;
$mode ||= 1;
if(!defined($k)){
	if($mode==1){$k=5}
	else{$k=15}
}
if(!defined($reads) || !defined($genome) || !defined($mat) || defined($help)){
	&usage;
	exit 0;
}

if(!defined($d)){
	$d="-d";
}else{
	$d="";
}
if(defined($c)){
	$c="-c";
}else{
	$c="";
}

my $pwd = `pwd`;
$pwd =~ s/\s+//g;

open OUT, ">$prj\_runnovel_miRNA.sh" or die $!;
print OUT "echo Start Time:\ndate\n";
if(!(-e "$pwd/$prj")){
	print OUT "mkdir -p $pwd/$prj\n";
}

if(!(-e "$pwd/$prj/ref")){
	print OUT "mkdir -p $pwd/$prj/ref\n";
}
print OUT "ln -sf ",get_ful_path($mat)," $pwd/$prj/ref/novel_mature_ref.fa\n";
my $genome_name = basename($genome);
print OUT "ln -sf ",get_ful_path($genome),"* $pwd/$prj/ref\n";
my $reads_name = basename($reads);
print OUT "ln -sf ",get_ful_path($reads)," $pwd/$prj/$reads_name\n\n";

print OUT "echo ================== Run predict novel miRNA ==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "export MIREVO=$miREvo\n";
print OUT "export PERL5LIB=$miREvo:\$PERL5LIB\n";
print OUT "export PATH=$srnatoolscli:\$PATH\n";
print OUT "export PATH=$bowtie:\$PATH\n";
print OUT "cd $pwd/$prj/\n";
print OUT "$Bin/pipline_predict.sh -o $prj.predict -i $reads_name -r ref/$genome_name -M ref/novel_mature_ref.fa -m $mode -k $k -p $threads $c -g $gear\n\n";

print OUT "echo ================== Run quantifier_gb novel miRNA ==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "$perlExec $Bin/novel_hp_mat_remove_dup.pl $prj.predict/predict.result.csv ref/predict_hairpin.fa ref/predict_mature.fa ref/predict_hairpin.pos\n";
#print OUT "awk \'{if(/^>/){split(\$1,a,\"_\");print \">$prj.novel_\"a[2]}else{gsub(\"T\",\"U\");gsub(\"t\",\"u\");print}}\' $prj.predict/predict.hairpin.fa >ref/predict_hairpin.fa\n";
#print OUT "awk \'{if(/^>/){split(\$1,a,\"_\");printf(\">$prj.novel_\"a[2]);if(/*\$/){printf(\"*\")};printf(\"\\n\")}else{gsub(\"T\",\"U\");gsub(\"t\",\"u\");print}}\' $prj.predict/predict.mature.fa >ref/predict_mature.fa\n";
#print OUT "awk \'{split(\$1,a,\"_\");print \"$prj.novel_\"a[2]\"\\t\"\$2}\' $prj.predict/predict.mirna.pos >ref/predict_hairpin.pos\n";
print OUT "export PATH=$mirdeep2:\$PATH\n";
print OUT "#export PATH=$mirdeep2/lib:\$PERL5LIB\n";
#print OUT "$perlExec $Bin/../bin/quantifier_gb.pl -p ref/predict_hairpin.fa -m ref/predict_mature.fa -r $reads_name -y $prj.novel -g $mismatch -T $threads $d\n";
print OUT "$perlExec $Bin/quantifier_gb_v2.pl -p ref/predict_hairpin.fa -m ref/predict_mature.fa -r $reads_name -y $prj.novel -g $mismatch -T $threads \n\n";

print OUT "echo ================== Run $prj.novel miRNA output ==================\n";
print OUT "export PATH=$srnatoolscli:$srnatoolscli/bin:$ViennaRNA2/bin:\$PATH\n";
print OUT "$perlExec $Bin/cvs2result.pl miRNAs_expressed_all_samples_$prj.novel.csv ref/predict_hairpin.fa ref/predict_mature.fa ref/predict_hairpin.pos $prj.novel\n";
print OUT "rm -rf $prj.novel/image/legend.txt\n";
print OUT "awk '{if(\$4>0){print \$1\"\\t\"\$3}}' miRNAs_expressed_all_samples_$prj.novel.csv >$prj.novel/hairpin_mature.pairs\n";
print OUT "head -1 $prj.novel/miRNAs_expressed_all_samples_$prj.novel.csv|awk \'{num=(NF-4)/2;printf(\"miRNA\");for(i=1;i<=num;i++){printf(\"\\t\"\$(i+4))};printf(\"\\n\");}\' >$prj.novel/mature.readcount\n";
print OUT "awk \'{if(NR>1 && \$1!~/*\$/){num=(NF-4)/2;printf(\$1\"\\t\"\$2);for(i=1;i<=num;i++){printf(\"\\t\"\$(i+4))};printf(\"\\n\")}}\' $prj.novel/miRNAs_expressed_all_samples_$prj.novel.csv|sort -k 1,1 -k 2nr,2|awk \'{if(\$1!=name){num=NF-2;printf(\$1);for(i=1;i<=num;i++){printf(\"\\t\"\$(i+2));}printf(\"\\n\");name=\$1}}\' >>$prj.novel/mature.readcount\n";
print OUT "head -1 $prj.novel/miRNAs_expressed_all_samples_$prj.novel.csv|awk \'{num=(NF-4)/2;printf(\"miRNA\");for(i=1;i<=num;i++){printf(\"\\t\"\$(i+4))};printf(\"\\n\");}\' >$prj.novel/star.readcount\n";
print OUT "awk \'{if(NR>1 && \$1~/*\$/){num=(NF-4)/2;printf(\$1\"\\t\"\$2);for(i=1;i<=num;i++){printf(\"\\t\"\$(i+4))};printf(\"\\n\")}}\' $prj.novel/miRNAs_expressed_all_samples_$prj.novel.csv|sort -k 1,1 -k 2nr,2|awk \'{if(\$1!=name){num=NF-2;printf(\$1);for(i=1;i<=num;i++){printf(\"\\t\"\$(i+2));}printf(\"\\n\");name=\$1}}\' >>$prj.novel/star.readcount\n";
print OUT "awk \'{if(/^>/){title=\$1}else{if(title){if(/total read count/ && (\$4>0)){print title;title=\"\";marker=1}else{marker=0}}if(marker){print}}}\' expression_analyses/expression_analyses_$prj.novel/miRBase.mrd >$prj.novel/miRBase.mrd\n\n";

print OUT "echo ================== Run Stat novel miRNA ==================\n";
print OUT "echo Start Time:\ndate\n";
print OUT "$perlExec $Bin/genebwt12count.pl -i expression_analyses/expression_analyses_$prj.novel/*_mapped.bwt.ka -r expression_analyses/expression_analyses_$prj.novel/precursor.converted -t novel_miRNA -o $prj.novel -u -s -W\n";
print OUT "awk \'{if(NR==1){print \"Types\\t\"\$0}else if(NR==2){print \"Mapped mature\\t\"\$0}}\' $prj.novel/novel_miRNA.mapmat.stat >$prj.novel/novel_miRNA.map.stat\n";
#edit by JC on 2015/12/23 同known miRNA============
print OUT "$perlExec $Bin/stat_known_miRNA_pre.pl $prj.novel/miRBase.mrd >$prj.novel/novel_miRNA.mapref.stat\n";
#=================================================
print OUT "awk \'{if(NR==2){print \"Mapped star\\t\"\$0}}\' $prj.novel/novel_miRNA.mapstar.stat >>$prj.novel/novel_miRNA.map.stat\n";
print OUT "awk \'{if(NR==2){print \"Mapped hairpin\\t\"\$0}}\' $prj.novel/novel_miRNA.mapref.stat >>$prj.novel/novel_miRNA.map.stat\n";
print OUT "awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=int(\$i+0.5)}printf(\"Mapped uniq sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"int(\$i+0.5))}printf(\"\\n\");}}\' $prj.novel/novel_miRNA.uc.stat >>$prj.novel/novel_miRNA.map.stat\n";
#print OUT "awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf(\"Mapped uniq sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"\$i)}printf(\"\\n\");}}\' $prj.novel/novel_miRNA.uc.stat >>$prj.novel/novel_miRNA.map.stat\n";
print OUT "awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=int(\$i+0.5)}printf(\"Mapped total sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"int(\$i+0.5))}printf(\"\\n\");}}\' $prj.novel/novel_miRNA.rc.stat >>$prj.novel/novel_miRNA.map.stat\n";

## 2014-2-10 add  for unifying "hairpin.fa" and last "*map.stat"
print OUT "cd $prj.novel\n";
#print OUT "$perlExec /PUBLIC/source/RNA/smallRNA/version3/bin/update_mapstat.pl hairpin.fa  novel_miRNA.map.stat  # for unify <hairpin.fa> <*.map.stat>\n";
print OUT "$RNAfold < hairpin.fa -noPS > hairpin.str\n";
print OUT "cd ../\n";

print OUT "$perlExec $Bin/bwt12collapse.pl expression_analyses/expression_analyses_$prj.novel/*_mapped.bwt.k1 >$prj.novel/novel_miRNA.map.fas\n";
print OUT "cp expression_analyses/expression_analyses_$prj.novel/*.unmap.fas $prj.novel/novel_miRNA.unmap.fas\n";
print OUT "$perlExec $Bin/bwt_base_bias_plot.pl -i expression_analyses/expression_analyses_$prj.novel/*_mapped.bwt.k1\n";
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

