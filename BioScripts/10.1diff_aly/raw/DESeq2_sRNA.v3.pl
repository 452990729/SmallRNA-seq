#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

sub usage
{
	print STDERR <<USAGE;
===============================================================================
Description     run DESeq2 and plotting
Version:1.0
2013-11-20      zhaohui\@novogene.cn

Options

-i  <s> : name of input file1,the reads matrix [merged.readcount]
-n1 <s> : treat group name
-n2 <s> : control group name
-a <s> : treat group samples, splitd by ":"
-b <s> : control group samples, splitd by ":"
-o  <s> : output_file,the outcome matrix, default='./DESeq2.out'
-op <s> : outputdir of pictures,default='./'
-p  <f> : padj-value threshold, default=0.05
-sm <s> : sharingMode of estimate Dispersions
-ft <s> : fit type of estimate Dispersions
-lf <s> : locfunc of estimate Sizefactors
===============================================================================
USAGE
}
my ($pv, $Fc, $sm, $ft, $lf, $in, $o, $Op, $n1, $n2, $group1, $group2, $help);
GetOptions(
        "h|?|help"=>\$help,
		"p=f"=>\$pv,
		"f=i"=>\$Fc,
		"sm=s"=>\$sm,
		"ft=s"=>\$ft,
		"lf=s"=>\$lf,
		"i=s"=>\$in,
		"o=s"=>\$o,
		"op=s"=>\$Op,
        "n1=s"=>\$n1,
        "n2=s"=>\$n2,
		"a=s"=>\$group1,
		"b=s"=>\$group2,
);
if(!defined($in) ||defined($help)){
        &usage;
        exit 0;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $R_v303 = $Config->{srnaenv}->{R_v303};

my $dir = `pwd`;
$dir =~ s/\n//;

$pv ||= 0.05;
$sm ||= 'maximum';
$ft ||= 'parametric';
$lf ||= 'median';
$o ||= "$dir/DESeq2.out";
$Op ||= "$dir";
mkdir $Op if (!-e $Op);
$n1 ||= 'group1';
$n2 ||= 'group2';
my $gname=$n1."vs".$n2;
my $readcount1=$n1."_readcount";
my $readcount2=$n2."_readcount";
my $R =<< "END";
#==============================================================================
#load data
#
library('reshape2')
library('DESeq2')
library("ggplot2")
library("pheatmap")
x=read.table("$in",head=T)
samplenames1<-unlist(strsplit("$group1",":"))
samplenames2<-unlist(strsplit("$group2",":"))
x=x[,c("sRNA",as.vector(samplenames1),as.vector(samplenames2))]
x<-x[order(x[,1]),]
A=length(samplenames1)
B=length(samplenames2)
group=factor(c(rep("$n1",times=A),rep("$n2",times=B)))
groups=c("$n1","$n2")
#------------------------------------------------------------------------------

## estimate parameter
sRNA_ID<-x[,1]
x_change<-x[,-1]
rownames(x_change)<-sRNA_ID
colData<-data.frame(condition=group)
row.names(colData)<-colnames(x_change)
dds<-DESeqDataSetFromMatrix(countData=x_change,colData=colData,design=~condition)
colData(dds)\$condition<-factor(colData(dds)\$condition,levels=groups)
dds <- DESeq(dds)
res <- results(dds)
res<-as.data.frame(res)
colnames(res)[5]<-"pval"
res\$log2FoldChange<--res\$log2FoldChange

####plot picture with DESeq2
##MA plot
png("$Op/$gname.MA.png",type="cairo-png")
plotMA(dds)
dev.off()
pdf("$Op/$gname.MA.pdf")
plotMA(dds)
dev.off()

##dispersion plot
png("$Op/$gname.dispersion.png",type="cairo-png")
plotDispEsts(dds)
dev.off()
pdf("$Op/$gname.dispersion.pdf")
plotDispEsts(dds)
dev.off()

#### calculate normalize parameter using DESeq
library(DESeq)
cds <- newCountDataSet(x_change, group)
cds <- estimateSizeFactors(cds) 
cds_norm<-counts(cds,T)
cds_norm<-as.data.frame(cds_norm)
cds_norm\$baseMeanA<-apply(cds_norm[,as.vector(samplenames1),drop=FALSE],1,mean)
cds_norm\$baseMeanB<-apply(cds_norm[,as.vector(samplenames2),drop=FALSE],1,mean)
Res<-cbind(cds_norm[,c(dim(cds_norm)[2]-1,dim(cds_norm)[2])],res[,c(2,5,6)])
Res\$significant<-Res\$padj<$pv

#------------------------------------------------------------------------------
#output files
sRNA_ID<-as.data.frame(sRNA_ID)
colnames(sRNA_ID)<-"sRNA"
out <- cbind(sRNA_ID,Res)
out[,c("log2FoldChange","pval","padj")]<-signif(out[,c("log2FoldChange","pval","padj")],5)
outSorted <- out[order(out\$padj),]
k=na.omit(outSorted)

out3 <- subset(outSorted,significant=='TRUE')
out_up<-subset(out3,out3\$log2FoldChange>0)
out_down<-subset(out3,out3\$log2FoldChange<0)
diff <- out3[,1]
diff_up<-out_up[,1]
diff_down<-out_down[,1]

diffsRNA <- out3[,c("sRNA","baseMeanA","baseMeanB","log2FoldChange","pval","padj")]
colnames(diffsRNA)<-c("sRNA","$readcount1","$readcount2","log2FoldChange","pval","padj")
diffsRNA_up<-out_up[,c("sRNA","baseMeanA","baseMeanB","log2FoldChange","pval","padj")]
colnames(diffsRNA_up)<-colnames(diffsRNA)<-c("sRNA","$readcount1","$readcount2","log2FoldChange","pval","padj")
diffsRNA_down<-out_down[,c("sRNA","baseMeanA","baseMeanB","log2FoldChange","pval","padj")]
colnames(diffsRNA_down)<-c("sRNA","$readcount1","$readcount2","log2FoldChange","pval","padj")
write.table(diffsRNA, file="$Op/$gname.DE.xls", sep="\t", quote=F, row.name=F)
write.table(diffsRNA_up, file="$Op/$gname.DE_up.xls", sep="\t", quote=F, row.name=F)
write.table(diffsRNA_down, file="$Op/$gname.DE_down.xls", sep="\t", quote=F, row.name=F)
write.table(diff, file="$Op/$gname.DElist.txt", sep="\t", quote=F, row.name=F, col.name=F)
write.table(diff_up, file="$Op/$gname.DElist_up.txt", sep="\t", quote=F, row.name=F, col.name=F)
write.table(diff_down, file="$Op/$gname.DElist_down.txt", sep="\t", quote=F, row.name=F, col.name=F)
out<-out[,c("sRNA","baseMeanA","baseMeanB","log2FoldChange","pval","padj","significant")]
colnames(out)<-c("sRNA","$readcount1","$readcount2","log2FoldChange","pval","padj","significant")
write.table(out, file="$Op/$gname.Differential_analysis_results.xls", sep="\t", quote=F, row.name=F)
#==============================================================================

## -------------------------------- Volcanno plot ------------------------------------------
setwd("$Op")
### k:all result without NA
### k_up,k_down,k_none

### code of zhaohui
# k_up<-subset(k,significant=="TRUE" & log2FoldChange>0)
# k_down<-subset(k,significant=="TRUE" & log2FoldChange<0)
# k_none<-subset(k,significant=="FALSE")
# up<-dim(k_up)[1]
# down<-dim(k_down)[1]
# total<-up+down
# uplable=paste("up:",up)
# downlable=paste("down:",down)
# k_up\$sig<-uplable
# k_down\$sig<-downlable
# k_none\$sig<-"FALSE"
# k<-rbind(k_up,k_down,k_none)
# logpadj<--log10(k\$padj)
# step<-floor(min(max(logpadj),300)/3)
# 
# p<-ggplot(k)+geom_point(aes(x=log2FoldChange,y=-log10(padj),color=sig),size=1.5)+ggtitle("$n1 vs $n2")+xlab(bquote(paste(log[2],"(fold change)",sep=""))) + ylab(bquote(paste(-log[10],"(padj)",sep="")))
# p<-p+scale_color_manual("",breaks=c(uplable,downlable,NA),values=c("green4","blue","red"))
# p<-p+geom_hline(yintercept=-log10($pv),linetype="dotdash",size=0.4)
# if (step>50) {
# p<-p+scale_y_continuous(breaks=c(-log10($pv),seq(step,min(max(logpadj),300),step)),labels=c(round(-log10($pv),1),seq(step,min(max(logpadj),300),step)))
# }
# if (step<=50) {
# p<-p+scale_y_continuous(breaks=c(-log10($pv),seq(0,min(max(logpadj),300),step)),labels=c(round(-log10($pv),1),seq(0,min(max(logpadj),300),step)))
# }
# ggsave(plot=p,file="$gname.Volcanoplot.pdf")
# ggsave(plot=p,file="$gname.Volcanoplot.png",type='cairo-png')

### code of me, 2014.1.1
## add up or down information
for(i in 1:nrow(k))
{
	if(k[i, "significant"]=="TRUE" & k[i, "log2FoldChange"]>0)
		k[i, "type"]="up"
	else if (k[i, "significant"]=="TRUE" & k[i, "log2FoldChange"]<0)
		k[i, "type"]="down"
	else
		k[i, "type"]="false"
}

upNum=nrow(subset(k, type=="up"))
downNum=nrow(subset(k, type=="down"))
falseNum=nrow(subset(k, type=="false"))

## function()
sampleVolcano <- function(dat,titles){
	p=ggplot(dat, aes(x=log2FoldChange, y=-log10(padj), colour=type)) + geom_point() + coord_cartesian(xlim=c(-5, 5)) + geom_hline(yintercept=-log10($pv), linetype=2) + labs(title=titles, x=bquote(paste(log[2], "(fold change)", sep="")), y=bquote(paste(-log[10],"(padj)",sep="")))

	## key, custom colour!
	p=p + scale_colour_manual(limits=c("up","down","false"), values=c("red", "green", "blue"), breaks=c("up", "down"), labels=c(paste("up:", upNum), paste("down:", downNum)))

	ggsave(filename="$gname.Volcanoplot.pdf", plot=p)	
	ggsave(filename="$gname.Volcanoplot.png", type="cairo-png", plot=p)
}

sampleVolcano(k, "$gname")

#==============================================================================
END
open FILE,">$Op/run_DESeq2.Rscript";
print FILE $R;
close FILE;
open R,"|$R_v303/R --vanilla --slave" or die $!;
print R $R;
close R;
