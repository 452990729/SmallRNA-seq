use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

sub usage{
print STDERR <<USAGE;
==================================================================================
example:
perl $0 --total S_1:../3map/prj/S_1.mapping.stat,S_2:../3map/prj/S_2.mapping.stat --unsense ../4k_miRNA/prj/prj.known/known_miRNA.uc.stat:../4k_miRNA/prj/prj.known/known_miRNA.rc.stat,../5ncRNA/prj/output/uc.stat:../5ncRNA/prj/output/rc.stat,../6repeat/repeat.uc.stat:../6repeat/repeat.rc.stat,../7NAT/prj/output/NAT.uc.stat:../7NAT/prj/output/NAT.rc.stat,../9n_miRNA/prj/prj.novel/novel_miRNA.uc.stat:../9n_miRNA/prj/prj.novel/novel_miRNA.rc.stat,../9_plus_plant_TAS/map/output/TAS.uc.stat:../9_plus_plant_TAS/map/output/TAS.rc.stat --sense ../8gene/prj/output/uc.stat:../8gene/prj/output/rc.stat
Usage: perl $0 [options]
options:
[mandatory parameters]
	--total sample1:mapping.stat,sample2:mapping.stat
		List of each sample's name and corresponeding mapping stat files, seperated by ","
	--unsense known_miRNA.uc.stat:known_miRNA.rc.stat,ncRNA/uc.stat:ncRNA/rc.stat,repeat.uc.stat:repeat.rc.stat,novel_miRNA.uc.stat:novel_miRNA.rc.stat
		List of "uc.stat:rc.stat" of each types(from known_miRNA to novel_miRNA, regardless of stranded type), seperated by ",",these types only display it's total mapped reads
[optional parameters]
	--sense	exon.uc.stat:exon.rc.stat,intron.uc.stat:intron.rc.stat
		List of "uc.stat:rc.stat" of each types(from known_miRNA to novel_miRNA, considering stranded type), seperated by ",",these types will display it's sense and antisense mapped reads
	-h|?|help
		Show this help
===============================================================================
USAGE
}

my ($totals,$unsenses,$senses,$help);
GetOptions(
	"h|?|help"=>\$help,
	"total=s"=>\$totals,
	"unsense=s"=>\$unsenses,
	"sense=s"=>\$senses,
);

if(defined($help)||!defined($totals)||!(defined($unsenses)||defined($senses))){
	&usage;
	exit 0;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../..//BioModule/Settings/config.ini");
my $R_v312 = $Config->{srnaenv}->{R_v312};

my @sample;
my %uthash;
my %rthash;
my @total=split(",",$totals);
foreach my $i(@total){
	my @tmp=split(":",$i);
	push(@sample,$tmp[0]);
	my $tmp1=`awk -F"\\t" '{if(NR==2){print \$2","\$4}}' $tmp[1]`;
	chomp($tmp1);
	my @tmp1=split(",",$tmp1);
	$rthash{$tmp[0]}=$tmp1[0];
	$uthash{$tmp[0]}=$tmp1[1];
#	print "$tmp[0]\t$rthash{$tmp[0]}\n";
#	print "$tmp[0]\t$uthash{$tmp[0]}\n";
}

my ($ut,$rt)=("c\(\"","c\(\"");
@sample=sort { $a cmp $b}  @sample;
for my $i(0..$#sample){
	if($i==$#sample){
		$ut.=$uthash{$sample[$i]}."\"\)";
		$rt.=$rthash{$sample[$i]}."\"\)";
	}else{
		$ut.=$uthash{$sample[$i]}."\",\"";
		$rt.=$rthash{$sample[$i]}."\",\"";
	}
}
#print "$ut\n$rt\n";

my $uc_category="sample.category_uc.txt";
my $rc_category="sample.category_rc.txt";
open(OUT1,">$uc_category");
open(OUT2,">$rc_category");
print OUT1 "Types\t",join("\t",@sample),"\n";
print OUT2 "Types\t",join("\t",@sample),"\n";

if(defined($unsenses)){
	my @unsense=split(",",$unsenses);
	foreach my $i(@unsense){
		my @tmp=split(":",$i);
		my $tmp0=`awk '{if(NR%3==2){print}}' $tmp[0]`;
		my $tmp1=`awk '{if(NR%3==2){print}}' $tmp[1]`;
		print OUT1 "$tmp0";
		print OUT2 "$tmp1";
	}
}

if(defined($senses)){
	my @sense=split(",",$senses);
	foreach my $i(@sense){
		my @tmp=split(":",$i);
		my $tmp0=`awk '{if(NR%3!=2 && NR!=1){print}}' $tmp[0]`;
		my $tmp1=`awk '{if(NR%3!=2 && NR!=1){print}}' $tmp[1]`;
		print OUT1 "$tmp0";
		print OUT2 "$tmp1";
	}

}

close(OUT1);
close(OUT2);

open(IN1,"$uc_category");
my $uc_category_full="category_uc_full.txt";
#my $uc_category_full=join("-",@sample).".category_uc_full.txt";
open(OUT1,">$uc_category_full");
my $header=<IN1>;
chomp($header);
my @heads=split("\t",$header);
print OUT1 $heads[0];
foreach my $i(1..$#heads){
	print OUT1 "\t$heads[$i]\t$heads[$i](percent)";
}
print OUT1 "\n";
print OUT1 "total";
for my $i(@sample){
	print OUT1 "\t",$uthash{$i},"\t100.00%";
}
print OUT1 "\n";
my %sum_u;
while(<IN1>){
	chomp;
	my @tmp=split("\t",$_);
	print OUT1 $tmp[0];
	foreach my $i(1..$#tmp){
		$sum_u{$sample[$i-1]}+=int($tmp[$i]+0.5);
		print OUT1 "\t",int($tmp[$i]+0.5),"\t",sprintf("%0.2f",$tmp[$i]*100/$uthash{$sample[$i-1]}),"%";
	}
	print OUT1 "\n";
}
print OUT1 "other";
for my $i(@sample){
	print OUT1 "\t",$uthash{$i}-$sum_u{$i},"\t",sprintf("%0.2f",($uthash{$i}-$sum_u{$i})*100/$uthash{$i}),"%";
}
print OUT1 "\n";
close(IN1);
close(OUT1);


open(IN2,"$rc_category");
my $rc_category_full="category_rc_full.txt";
#my $rc_category_full=join("-",@sample).".category_rc_full.txt";
open(OUT2,">$rc_category_full");
$header=<IN2>;
chomp($header);
@heads=split("\t",$header);
print OUT2 $heads[0];
foreach my $i(1..$#heads){
        print OUT2 "\t$heads[$i]\t$heads[$i](percent)";
}
print OUT2 "\n";
print OUT2 "total";
for my $i(@sample){
	print OUT2 "\t",$rthash{$i},"\t100.00%";
}
print OUT2 "\n";
my %sum_r;
while(<IN2>){
        chomp;
        my @tmp=split("\t",$_);
        print OUT2 $tmp[0];
        foreach my $i(1..$#tmp){
                $sum_r{$sample[$i-1]}+=int($tmp[$i]+0.5);
                print OUT2 "\t",int($tmp[$i]+0.5),"\t",sprintf("%0.2f",$tmp[$i]*100/$rthash{$sample[$i-1]}),"%";
        }
        print OUT2 "\n";
}
print OUT2 "other";
for my $i(@sample){
        print OUT2 "\t",$rthash{$i}-$sum_r{$i},"\t",sprintf("%0.2f",($rthash{$i}-$sum_r{$i})*100/$rthash{$i}),"%";
}
print OUT2 "\n";
close(IN2);
close(OUT2);

`rm $uc_category $rc_category`;
my $pwd=`pwd`; $pwd=~s/\n//g;
my $R =<< "END";
#==============================================================================
library('ggplot2')
setwd("$pwd");
ut<-$ut
rt<-$rt
uc<-read.delim("$uc_category_full",sep="\\t",header=TRUE)
n_list<-c(1:length(colnames(uc)))
n_list<-c(1,n_list[n_list%%2==0])
uc<-uc[c(2:length(rownames(uc))),c(n_list)]
n<-length(colnames(uc))
for(i in 2:n){
	pielabels<-c()
        pieval<-uc[,i]
        for(j in 1:length(rownames(uc))){
                pielabels<-c(pielabels,paste(uc[j,1],"(",uc[j,i],", ",round(as.numeric(uc[j,i])*100/as.numeric(ut[i-1]),2),"%)",sep=""))
        }
        pieval[pieval==0]<-NA
        dat=data.frame(pieval,pielabels)
        p<-ggplot(na.omit(dat),aes(x='',y=pieval,fill=pielabels))+geom_bar(width = 1)+scale_y_continuous(breaks=c())+coord_polar(theta = "y")+theme(legend.title=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank())+ggtitle(paste("Annotation of Uniq reads\\n","(", colnames(uc)[i],")",sep=""))+scale_fill_hue(l=40)
        ggsave(filename=paste(colnames(uc)[i],".category_uc_pie.png",sep=""),type='cairo',dpi=300)
        ggsave(filename=paste(colnames(uc)[i],".category_uc_pie.pdf",sep=""))
	dev.off()
}
rc<-read.delim("$rc_category_full",sep="\\t",header=TRUE)
n_list<-c(1:length(colnames(rc)))
n_list<-c(1,n_list[n_list%%2==0])
rc<-rc[c(2:length(rownames(rc))),c(n_list)]
n<-length(colnames(rc))
for(i in 2:n){
	pielabels<-c()
	pieval<-rc[,i]
	for(j in 1:length(rownames(rc))){
		pielabels<-c(pielabels,paste(rc[j,1],"(",rc[j,i],", ",round(as.numeric(rc[j,i])*100/as.numeric(rt[i-1]),2),"%)",sep=""))
	}
	pieval[pieval==0]<-NA
	dat=data.frame(pieval,pielabels)
	p<-ggplot(na.omit(dat),aes(x='',y=pieval,fill=pielabels))+geom_bar(width = 1)+scale_y_continuous(breaks=c())+coord_polar(theta = "y")+theme(legend.title=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background = element_blank())+ggtitle(paste("Annotation of Total reads\\n","(", colnames(rc)[i],")",sep=""))+scale_fill_hue(l=40)
	ggsave(filename=paste(colnames(rc)[i],".category_rc_pie.png",sep=""),type='cairo',dpi=300)
	ggsave(filename=paste(colnames(rc)[i],".category_rc_pie.pdf",sep=""))
	dev.off()
}

#==============================================================================
END

open FILE, ">graph.R";
print FILE $R;
close FILE;
system "$R_v312/Rscript graph.R";

