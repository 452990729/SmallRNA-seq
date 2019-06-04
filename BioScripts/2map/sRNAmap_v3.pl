#2013.07.24.by zhangyu
#2013.07.26.by wangshaobin 
#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
#my $Bin="/PUBLIC/source/RNA/smallRNA/version3/2map";
use Config::Tiny;

sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description     Creat sh to run sRNA mapping analysis.
Version:1.0  
Usage: perl $0 [options]
       qsub -cwd -V -l vf=1g -l p=5 runsRNA_mapping.sh
Options
[mandatory parameters]
                -r  <s>     :  Reference genome or transcript sequence
                -q  <s>     :  QC directory
                -p  <s>     :  List of sample name(3 character) corresponding to reads, seperated by ","
                -o  <s>     :  Output directory
                -h|?|help   :  Show this help
[optional parameters]
				-n  [int]   :  The number of chrosome to be shown, use for genome
                -v  [int]   :  Number of allowed mismatches when mapping reads to precursors, default 0
                -T  [int]   :  Number of alignment threads to launch during bowtie mapping, default 5
                -w  [int]   :  The window size of each point, default 1000
                -s  [int]   :  The step size, default 1000
=========================================================================
USAGE
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $python = $Config->{srnaenv}->{python_v276};
my $bowtie = $Config->{software}->{bowtie1};
my $samtools = $Config->{software}->{samtools};
my $circos = $Config->{software}->{circos};
my $R_v2153 = $Config->{srnaenv}->{R_v2153};
my $LIBRARY_PATH = $Config->{srnaenv}->{LIBRARY_PATH};

my ($help,$refer,$qc,$samples,$out,$number,$mismatchs,$threads,$window,$step);
GetOptions(
        "h|?|help"=>\$help,
        "r=s"=>\$refer,
        "q=s"=>\$qc,
        "p=s"=>\$samples,
        "o=s"=>\$out,
        "n=i"=>\$number,
        "v=i"=>\$mismatchs,
        "T=i"=>\$threads,
        "w=i"=>\$window,
        "s=i"=>\$step,

);

$mismatchs = 0 if(!defined($mismatchs));
$threads ||= 5;
$window ||= 1000;
$step ||= 1000;

if(!defined($refer) || !defined($qc) || !defined($samples) || !defined($out) ||defined($help)){
        &usage;
        exit 0;
}

my $pwd = `pwd`;
$pwd =~ s/\s+//g;

open OUT, ">runsRNA_mapping.sh" or die $!;

if(!-e "$pwd/$out")
{
	mkdir "$pwd/$out";
}
if(!-e "$pwd/$out/tmp")
{
	mkdir "$pwd/$out/tmp";
}

#$refer=basename($refer);
print OUT "echo Start Time:\ndate\n";
my @sample=(split /,/, $samples);
my @read;
foreach my $i(@sample){
        my $j="$qc/$i/clean_data/$i\_remain_uniq.fa";
        push @read,$j;
}
#mapping
print OUT "$perlExec -e  \'print \"Sample\\tTotal sRNA\\tMapped sRNA\\t\\\"+\\\" Mapped sRNA\\t\\\"-\\\" Mapped sRNA\\n\"\' >$pwd/$out/reference.mapping.stat\n";
foreach my $i(0..$#sample)
{
		open Each, ">$pwd/$out/$sample[$i]_map.sh";
		if(-e "$read[$i]")
		{
			print Each "echo \"===========================mapping for $sample[$i]:========================\"\n";
			#print Each "ln -s ",get_ful_path($read[$i])," $pwd/$out/tmp/\n";
			print Each "awk \'{if(/^>/){split(substr(\$1,2),a,\"(\");split(a[2],b,\")\");id=\"$sample[$i]_\"i++;for(j=1;j<=b[1];j++){print \">\"id\"_x\"b[1]\"_\"j\"\\n\"a[1]}}}\' $read[$i] >$pwd/$out/tmp/$sample[$i].total.fa\n";
			print Each "awk \'{if(/^>/){split(substr(\$1,2),a,\"(\");split(a[2],b,\")\");id=\"$sample[$i]_\"i++;print \">\"id\"_x\"b[1]\"\\n\"a[1]}}\' $read[$i] >$pwd/$out/tmp/$sample[$i].collapse.fa\n";
			print Each "$bowtie/bowtie -p $threads -v $mismatchs -k 1 $refer -f $pwd/$out/tmp/$sample[$i].total.fa --sam $pwd/$out/tmp/$sample[$i].sam 2>$pwd/$out/tmp/$sample[$i].mapping.log\n";
			print Each "$bowtie/bowtie -p $threads -v $mismatchs -k 1 $refer -f $pwd/$out/tmp/$sample[$i].total.fa $pwd/$out/tmp/$sample[$i].bwt --un $pwd/$out/tmp/$sample[$i].unmap.total.fa\n\n";

			print Each "echo \"=======================Stat mapping for $sample[$i]:=======================\"\n";
			print Each "$perlExec $Bin/bwt2sta-uniq.pl $pwd/$out/tmp/$sample[$i].bwt $pwd/$out/$sample[$i].mapping.stat $pwd/$out/tmp/$sample[$i].map.collapse.fa\n";
			print Each "awk \'{if(/^>/){split(\$1,a,\"_x\");uniq++;count=a[2];total+=count}else{uniqbase+=length(\$1);totalbase+=count*length(\$1)}}END{print \"Total small RNA\\t\"total\"\\t\"totalbase\"\\t\"uniq\"\\t\"uniqbase}\' $pwd/$out/tmp/$sample[$i].collapse.fa >>$pwd/$out/$sample[$i].mapping.stat\n";
			print Each "awk -F\"\\t\" \'BEGIN{sample=\"$sample[$i]\"}{if(/^Total small RNA/){total=\$2}else if(/^Total Mapped small RNA/){map=\$2}else if(/^Total Sense Mapped small RNA/){sense=\$2}else if(/^Total Antisense Mapped small RNA/){antisense=\$2}}END{map_percent=sprintf(\"%.2f%\",map/total*100);sense_percent=sprintf(\"%.2f%\",sense/total*100);antisense_percent=sprintf(\"%.2f%\",antisense/total*100);print sample\"\\t\"total\" (100.00%)\\t\"map\" (\"map_percent\")\\t\"sense\" (\"sense_percent\")\\t\"antisense\" (\"antisense_percent\")\"}' $pwd/$out/$sample[$i].mapping.stat >>$pwd/$out/reference.mapping.stat\n\n";

#distribution 
			if(defined($number))
			{
				print Each "echo \"=======================Genome distribution for $sample[$i]:=================\"\n";
				print Each "export PATH=$circos/bin:\$PATH\n";
                                print Each "export PATH=$R_v2153:\$PATH\n";
                                print Each "export LD_LIBRARY_PATH=$LIBRARY_PATH/R_v2153:$LIBRARY_PATH/HTSeq_Atlas:\$LD_LIBRARY_PATH\n";
				my $dir="$pwd/$out/tmp/$sample[$i]";
				print Each "if [ ! -e $dir ]; then mkdir -p $dir; fi \n";
				print Each "cd $pwd/$out/tmp/$sample[$i]\n";
				print Each "$samtools view -bS $pwd/$out/tmp/$sample[$i].sam  >$sample[$i].bam\n";
				print Each "$python $Bin/circos_density_v2.3.py --bam $sample[$i].bam --fa $refer --n $number --iv 10000 --r 10000 --output $sample[$i].circos.png \n";
				print Each "cp $pwd/$out/tmp/$sample[$i]/$sample[$i]*.png $pwd/$out/\n";
				print Each "cp $pwd/$out/tmp/$sample[$i]/$sample[$i]*.svg $pwd/$out/\n\n\n";
			}
		}
	close Each;
	print OUT "qsub -V -cwd -l vf=2G  $pwd/$out/$sample[$i]_map.sh","\n";
}

#print OUT "cat $pwd/$out/tmp/*.map.collapse.fa >$pwd/$out/all.map.collapse.fa\n";
close OUT;
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
