#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

sub usage
{
        print STDERR <<USAGE;
==================================================================================
Description  run DEGseq and plotting
2012-12-06  jiangxiaoxue\@novogene.cn

Options
<Required>
	-i <s>  :	meanscount.txt
	-r <s>	:	meanstpm.txt
	-a <s> : group1_name
	-b <s> : group2_name
<Optional>
	-o <s>  : 	 the output directory, default ./DEGseq.out
	-p <f>  : 	 the threshold of adjusted pvalue for diff analysis, default=0.01
	-f <f>  :  	 the threshold of |log2(Foldchange)| for diff analysis, default=2
===============================================================================
USAGE
}

my($rc, $tpm, $dir, $pval, $Fc, $group1, $group2, $help);
GetOptions(
        "h|?|help"=>\$help,
		"i=s"=>\$rc,
		"r=s"=>\$tpm,
		"o=s"=>\$dir,
		"p=f"=>\$pval,
		"f=f"=>\$Fc,
		"a=s"=>\$group1,
		"b=s"=>\$group2,
);


if((!defined($rc)) || (!defined($tpm))|| (!defined($group1))||(!defined($group2))|| (defined($help))){
        &usage;
        exit 0;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $R_v312 = $Config->{srnaenv}->{R_v312};

$pval ||= 0.01;
$Fc ||= 2;
my $pwd = `pwd`;
$pwd =~ s/\s+//;
$dir ||= "$pwd/DEGseq.out";

if(!(-e $dir)){
	`mkdir $dir`;
}

my $R =<< "END";
#==============================================================================
library("DEGseq")
library("edgeR")
library('reshape2')
library("ggplot2")
setwd("$dir")

############# normalization of readcount ##########################
Dat<- read.table("$rc", header=T, row.names=1)
Dat=Dat[,c("$group1","$group2")]
Dat<- Dat[order(rownames(Dat)),]
names<-colnames(Dat);
n<-length(colnames(Dat));
d<- DGEList(counts = Dat, group = names);
d.norm<- calcNormFactors(d);
write.table(d.norm\$samples,file="lib.size_norm.factors",sep = "\\t", quote=F);

Dat.nor<- Dat;
for(i in 1:n){
	Dat.nor[,i]<- Dat[,i]*1000000/((d.norm\$samples[i,2])*(d.norm\$samples[i,3]))
};
write.table(Dat.nor, file="$dir/all.counts.genes.matrix.nor", sep="\\t",quote=FALSE);

################## run diff_analysis ###########################
temp<-data.frame(GeneID=c(rownames(Dat.nor)), Chr=rep( "chr1", times=dim(Dat.nor)[1]), GeneStart=rep(1, times=dim(Dat.nor)[1]), GeneEnd=rep(1000, times=dim(Dat.nor)[1]), Status=rep("KNOWN", times=dim(Dat.nor)[1]), ExternalID=c(rownames(Dat.nor)) );
rownames(temp)<-rownames(Dat.nor);
countTable<-cbind(temp, Dat.nor);
write.table(countTable, file="$dir/countTable", quote=F, row.names=F, sep="\\t");
geneExpFile<-c("countTable")
geneExpMatrix1 <- readGeneExp(file = geneExpFile, geneCol = 1, valCol = c(6 + 1));
geneExpMatrix2 <- readGeneExp(file = geneExpFile, geneCol = 1, valCol = c(6 + 2));

DEGexp(geneExpMatrix1 = geneExpMatrix1, geneCol1 = 1, expCol1 = c(2), groupLabel1 = "names[array[i,1]]", geneExpMatrix2 = geneExpMatrix2, geneCol2 = 1, expCol2 = c(2), groupLabel2 = "names[array[i,2]]", method = "MARS", qValue=$pval , foldChange=$Fc , thresholdKind=5, outputDir="$dir",  normalMethod="none")

##################### Calculate TPM #######################################
tpm=read.table("$tpm",head=T)
tpm<-tpm[order(tpm[,1]),]
matrixTPM<-tpm[,c("$group1","$group2")]
colnames(matrixTPM)=paste(colnames(matrixTPM),"tpm",sep=".")
rc_tpm<-cbind(tpm[,1],Dat,matrixTPM)
colnames(rc_tpm)[c(1,2,3)]<-c("sRNA",colnames(Dat))
#colnames(rc_tpm)[c(1,2,3)]<-c("sRNA",paste(colnames(Dat),"readcount",sep="."))

### to extract diff list ##########################################
res<-read.table("output_score.txt", header=TRUE);
#res[,c("log2.Fold_change.","log2.Fold_change..normalized")]<- -res[,c("log2.Fold_change.","log2.Fold_change..normalized")]
res<-res[order(res[,1]),];
a=res[,4]>log2($Fc)|res[,4]<(-log2($Fc))
b=res[,9]<$pval
logic=a\&b
res[,10]=as.matrix(logic)
colnames(res)[10]<-"Signature"

out<-cbind(rc_tpm,res[,-1]);
colnames(out)[c(2,3,6,7)]<-c("$group1.readcount","$group2.readcount","$group1","$group2")
out[,c("log2.Fold_change.","p.value","q.value.Storey.et.al..2003.")]<-signif(out[,c("log2.Fold_change.","p.value","q.value.Storey.et.al..2003.")],5)
out <- out[order(out[,c("q.value.Storey.et.al..2003.")]),]

diff<-subset(out, out[14]==TRUE);
diff_up<-subset(diff,diff[,c("log2.Fold_change.")]>0);
diff_down<-subset(diff,diff[,c("log2.Fold_change.")]<0);
out_info<-cbind(out[,c(1,6,7)],out[,c("log2.Fold_change.","p.value","q.value.Storey.et.al..2003.","Signature")]);
write.table(out_info, file=paste("$group1","vs","$group2", ".Differential_analysis_results.xls", sep=""), quote=F, row.names=F, sep="\\t");

diff_info<-cbind(diff[,c(1,6,7)],diff[,c("log2.Fold_change.","p.value","q.value.Storey.et.al..2003.")]);
diff_up_info<-cbind(diff_up[,c(1,6,7)],diff_up[,c("log2.Fold_change.","p.value","q.value.Storey.et.al..2003.")]);
diff_down_info<-cbind(diff_down[,c(1,6,7)],diff_down[,c("log2.Fold_change.","p.value","q.value.Storey.et.al..2003.")]);
write.table(diff_info,file=paste("$group1","vs","$group2", ".DE.xls", sep=""), quote=F, row.names=F, sep="\\t");
write.table(diff_up_info,file=paste("$group1","vs","$group2", ".DE_up.xls", sep=""), quote=F, row.names=F, sep="\\t");
write.table(diff_down_info,file=paste("$group1","vs","$group2", ".DE_down.xls", sep=""), quote=F, row.names=F, sep="\\t");

diff_id<-diff[,1];
diff_up_id<-diff_up[,1];
diff_down_id<-diff_down[,1];
write.table(diff_id,file=paste("$group1","vs","$group2", ".DElist.txt", sep=""), quote=F, row.names=F, col.name=F, sep="\\t");
write.table(diff_up_id,file=paste("$group1","vs","$group2", ".DElist_up.txt", sep=""), quote=F, row.names=F, col.name=F, sep="\\t");
write.table(diff_down_id,file=paste("$group1","vs","$group2", ".DElist_down.txt", sep=""), quote=F, row.names=F, col.name=F, sep="\\t");

## ------------------------------------  volcano plot --------------------------------------------#
k<-out[,c(1,8,13,14)]
colnames(k)<-c("sRNA", "logFC","qval","signature")
k<-na.omit(k)

## add up or down information
for(i in 1:nrow(k))
{
	if(k[i, "signature"]=="TRUE" & k[i, "logFC"]>0)
		k[i, "type"]="up"
	else if (k[i, "signature"]=="TRUE" & k[i, "logFC"]<0)
		k[i, "type"]="down"
	else
		k[i, "type"]="false"
}

upNum=nrow(subset(k, type=="up"))
downNum=nrow(subset(k, type=="down"))
falseNum=nrow(subset(k, type=="false"))

## function()
sampleVolcano <- function(dat,titles){
	p=ggplot(dat, aes(x=logFC, y=-log10(qval), colour=type))+theme(panel.background=element_blank())+theme(axis.line=element_line(colour="black"))+theme(legend.position=c(0.9,0.85))+theme(panel.grid.major=element_line(colour=NA)) +theme(panel.grid.minor=element_line(colour=NA)) + theme(axis.text.x = element_text(face="bold",colour="black")) +theme(axis.text.y = element_text(face="bold",colour="black"))+ geom_point() + coord_cartesian(xlim=c(-5, 5)) + geom_hline(yintercept=2, linetype=2) + geom_vline(xintercept=c(-log($Fc,2),log($Fc,2)), linetype=2) + labs(title=titles, x=bquote(paste(log[2], "(fold change)", sep="")), y=bquote(paste(-log[10],"(qvalue)",sep="")))

	## key, custom colour!
	p=p + scale_colour_manual(limits=c("up","down","false"), values=c("red", "green", "blue"), breaks=c("up", "down"), labels=c(paste("up:", upNum), paste("down:", downNum))) + theme(legend.text=element_text(size=15)) + theme(legend.title=element_text(size=15))

	ggsave(filename=paste("$group1","vs","$group2", ".Volcanoplot.pdf", sep=""), plot=p)	
	ggsave(filename=paste("$group1","vs","$group2", ".Volcanoplot.png", sep=""), type="cairo-png", plot=p,dpi=300)
}

sampleVolcano(k, titles=paste("$group1"," vs ", "$group2"))
#==============================================================================
END

open(OUT,">$dir/DEGseq.Rscript");
print OUT $R;
close(OUT);

open R,"|$R_v312/R --vanilla --slave" or die $!;
print R $R;
close R;
