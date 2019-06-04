use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

my $usage=<<END;
Program:
 This program only map reads to human exon and intron;
History:
 2015/6/30 first release by jc
Eg:
 perl $0
	-fa repeat.unmap.fas
	-spe species three abbreviation in miRbase
	-dir output dir
	-species name of genome 
END

onfig = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $bowtie1 = $Config->{software}->{bowtie1};



my ($fa,$abbr,$dir,$species);
GetOptions(
	"fa=s" =>\$fa,
	"spe=s" =>\$abbr,
	"dir=s" =>\$dir,
	"species=s" =>\$species,
);

die $usage if (!$fa or !$abbr or !$dir);

open OUT,">>${abbr}_gene_map.sh";
my $exon_intron_dir = $Config->{database}->{hsaRef};
my $exon="$exon_intron_dir/$species/exon/hsa_${species}_exon.fa";
my $intron_1="$exon_intron_dir/$species/intron/hsa_${species}_intron_1.fa";
my $intron_2="$exon_intron_dir/$species/intron/hsa_${species}_intron_2.fa";
my $intron_3="$exon_intron_dir/$species/intron/hsa_${species}_intron_3.fa";
my $intron_4="$exon_intron_dir/$species/intron/hsa_${species}_intron_4.fa";
my $intron_5="$exon_intron_dir/$species/intron/hsa_${species}_intron_5.fa";
my $intron="$exon_intron_dir/$species/intron/hsa_${species}_intron.fa";;
my $run=<<END;
#
echo ================== Run exon mapping ==================
$bowtie/bowtie -p 10 -v 0 -k 1 -f $exon $fa $dir/exon.bwt --un $dir/exon.unmap.fas 2>$dir/exon.mapping.log
echo ================== Stat exon mapping  ==================
if [[ ! -z $dir/exon.bwt ]]; then
$perlExec $Bin/bwt12collapse.pl $dir/exon.bwt >$dir/exon.map.fas
cat $dir/exon.map.fas >$dir/gene.map.fas
$perlExec $Bin/genebwt12count.pl -i $dir/exon.bwt -r $exon -t exon -o $dir -u -s 
awk '{if(NR==1){print "Types\\t"\$0}else if(NR==2){print "Mapped reference\\t"\$0}}' $dir/exon.mapref.stat >$dir/exon.map.stat
awk '{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf("Mapped uniq sRNA\\t"total);for(i=2;i<=NF;i++){printf("\\t"\$i)}printf("\\n");}}' $dir/exon.uc.stat >>$dir/exon.map.stat
awk '{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf("Mapped total sRNA\\t"total);for(i=2;i<=NF;i++){printf("\\t"\$i)}printf("\\n");}}' $dir/exon.rc.stat >>$dir/exon.map.stat
fi
echo ================== Run intron mapping ==================
$bowtie/bowtie -p 10 -v 0 -k 1 -f $intron_1 $dir/exon.unmap.fas $dir/intron1.bwt --un $dir/intron1.unmap.fas 2>$dir/intron1.mapping.log
$bowtie/bowtie -p 10 -v 0 -k 1 -f $intron_2 $dir/intron1.unmap.fas $dir/intron2.bwt --un $dir/intron2.unmap.fas 2>$dir/intron2.mapping.log
$bowtie/bowtie -p 10 -v 0 -k 1 -f $intron_3 $dir/intron2.unmap.fas $dir/intron3.bwt --un $dir/intron3.unmap.fas 2>$dir/intron3.mapping.log
$bowtie/bowtie -p 10 -v 0 -k 1 -f $intron_4 $dir/intron3.unmap.fas $dir/intron4.bwt --un $dir/intron4.unmap.fas 2>$dir/intron4.mapping.log
$bowtie/bowtie -p 10 -v 0 -k 1 -f $intron_5 $dir/intron4.unmap.fas $dir/intron5.bwt --un $dir/intron5.unmap.fas 2>$dir/intron5.mapping.log
cat $dir/intron*.bwt >$dir/intron.bwt
echo ================== Stat intron mapping  ==================
if [[ ! -z $dir/intron.bwt ]]; then
$perlExec $Bin/bwt12collapse.pl $dir/intron.bwt >$dir/intron.map.fas
cat $dir/intron.map.fas >>$dir/gene.map.fas
$perlExec $Bin/genebwt12count.pl -i $dir/intron.bwt -r $intron -t intron -o $dir -u -s
awk '{if(NR==1){print "Types\\t"\$0}else if(NR==2){print "Mapped reference\\t"\$0}}' $dir/intron.mapref.stat >$dir/intron.map.stat
awk '{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf("Mapped uniq sRNA\\t"total);for(i=2;i<=NF;i++){printf("\\t"\$i)}printf("\\n");}}' $dir/intron.uc.stat >>$dir/intron.map.stat
awk '{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf("Mapped total sRNA\\t"total);for(i=2;i<=NF;i++){printf("\\t"\$i)}printf("\\n");}}' $dir/intron.rc.stat >>$dir/intron.map.stat
fi
$perlExec $Bin/paste_col.pl $dir/exon.uc.stat $dir/intron.uc.stat >$dir/uc.stat
$perlExec $Bin/paste_col.pl $dir/exon.rc.stat $dir/intron.rc.stat >$dir/rc.stat
echo End Time:
cp $dir/intron5.unmap.fas $dir/gene.unmap.fas
END
print OUT $run;
close OUT;
