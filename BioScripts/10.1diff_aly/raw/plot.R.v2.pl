#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);
use Config::Tiny;

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $R_v2153 = $Config->{srnaenv}->{R_v2153};

# Used for Plot density box, cluster 

my ($rc, $tpm, $mean_rc, $mean_tpm, $compare, $output, $cluster);
GetOptions(
		"rc=s"=>\$rc,
		"tpm=s"=>\$tpm,
		"compare=s"=>\$compare,
		"output=s"=>\$output,
		"cluster=s"=>\$cluster,
);
my $usage=<<END;
------------------------------------------------------------------------------
perl $0
-rc		readcount file
-tpm		tpm file
-compare	compare.txt
-output		outdir
-cluster	heat_sampel order
------------------------------------------------------------------------------
END
die $usage if (!$rc or !$tpm or !$compare or !$output);
$rc=abs_path($rc);
$tpm=abs_path($tpm);
$compare=abs_path($compare);
mkdir $output if(!-e $output);
$output=abs_path($output);

my $R =<<END;
# ==============================================================================
#load data
library("Vennerable")
library('reshape2')
library("ggplot2")
library("pheatmap")
setwd("$output")

# ==============================================================================
# load readcount
matrix_rc<-read.table("$rc",head=T)
matrix_rc<-matrix_rc[order(matrix_rc[,1]),]
colnames(matrix_rc)<-paste(colnames(matrix_rc),"readcount",sep=".")

# ==============================================================================
# load tpm
matrix_tpm<-read.table("$tpm",head=T)
matrix_tpm<-matrix_tpm[order(matrix_tpm[,1]),]
matrix_tpm<-matrix_tpm[-1]
colnames(matrix_tpm)<-paste(colnames(matrix_tpm),"tpm",sep=".")

# ==============================================================================
# merge tpm readcount
rc_tpm<- cbind(matrix_rc,matrix_tpm)
write.table(rc_tpm, file="Readcount_TPM.xls", sep="\t", quote=F, row.names=F)

# ==============================================================================
# tpm Interval
sRNA_num=length(matrix_tpm[1,])
TPM<-array(0:0,c(6,length(matrix_tpm[1,])))
TPM_number=rbind(colSums(matrix_tpm<=0.1), colSums(matrix_tpm<=0.3 & matrix_tpm>0.1), colSums(matrix_tpm<=3.57 & matrix_tpm>0.3),  colSums(matrix_tpm>3.57 & matrix_tpm<=15),  colSums(matrix_tpm<=60 & matrix_tpm>15),  colSums(matrix_tpm>60))
TPM_p=100*TPM_number/dim(matrix_tpm)[1]
for(i in 1:6)
{
	for(j in 1:sRNA_num)
	{
		TPM_p[i,j]=sprintf("%.2f",as.numeric(TPM_p[i,j]))
		TPM[i,j]<-paste(TPM_number[i,j],"(",TPM_p[i,j],"%)",sep="")
	}
}
TPM<-data.frame(TPM)
colnames(TPM)<-colnames(matrix_tpm)
TPM\$Interval<-c("0-0.1","0.1-0.3","0.3-3.57","3.57-15","15-60",">60")
TPM<-TPM[,c(sRNA_num+1,1:sRNA_num)]
names(TPM)[1]<-"TPM Interval"
out4<-TPM
write.table(out4,file="TPM_interval.xls",sep="\\t", quote=F, row.name=F)

# ==============================================================================
# Density and Box plot 
rp<-read.table("$tpm",head=T)
df<-log10(rp[-1]+1)
df<-melt(df)
colnames(df)<-c("Group","value")

## Density plot
p<-ggplot(df, aes(x=value, colour=Group, group=Group, fill=Group)) +  geom_density(alpha=0.2) + xlab("log10(TPM+1)") + ylab("Density") + labs(title="TPM density distribution")
pdf("TPM_density_distribution.pdf")
p
dev.off()
png("TPM_density_distribution.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
p
dev.off()

## Box plot
p<-ggplot(df, aes(Group, value)) + geom_boxplot(aes(fill = Group)) + xlab("") + ylab("log10(TPM+1)") + labs(title="TPM distribution")
pdf("TPM_boxplot.pdf")
p
dev.off()
png("TPM_boxplot.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
p
dev.off()

# ==============================================================================
# correlation plot
source("$Bin/corr_plot.r")
corr_plot("$tpm")

# ==============================================================================
# tpm hcluster heatmap: log10(tpm+1) 
compare<-read.table("$compare",head=F)
array<-compare[,1]
array<-as.vector(array)
IDdir<-paste(array,"/",array,".DElist.txt",sep="")   # The Place where error occur often because of no DE gene! 
dir<-paste(array,"/",array,".Differential_analysis_results.xls",sep="")   

un<-as.matrix("sRNA")
k<-as.matrix(rp[,1])
colnames(k)<-"sRNA"

for(i in 1:length(array))
{
	#a<-read.table(IDdir[i],head=F)  # The Place where error occur often because of no DE gene! 
	sRNA<-read.table(dir[i],header=T)  # The Place where error occur often because of no DE gene! 
	n=length(sRNA[1,])
	sRNA=sRNA[order(sRNA[,1]),]

	sig=as.matrix(sRNA[,n])
	colnames(sig)=array[i]
	k=cbind(k, sig)

	#if(is.na(file.info(IDdir[i])\$size)){
	if(file.info(IDdir[i])\$size==0){
		next
	}else{
		a=read.table(IDdir[i], header=F)
		un=rbind(un, as.matrix(a))
	}
	#un=as.matrix(un[-1,])
	#index=duplicated(un[,1])
	#union1=as.matrix(un[!index,])
	#union=subset(rp, sRNA %in% union1[,1])
}
un<-as.matrix(un[-1,])
index<-duplicated(un[,1])
union1<-as.matrix(un[!index,])
union_all<-subset(rp, sRNA %in% union1[,1])
#add0start
clusters = unlist(strsplit("sRNA\,$cluster",","))
union<-union_all[,clusters]
#add0end
write.table(union, file="DE_union_for_cluster", row.names=F, sep="\\t", quote=F, col.names=T)
write.table(k, file="$output/Signature.xls", row.names=F, sep="\\t", quote=F, col.names=T)

# ploting heatmap
rownames(union)<-union[,1]
union<-union[,-1]
union<-log10(union + 1)
#add start
#clusters = unlist(strsplit("$cluster",","))
#picname = paste(clusters,collapse="_")
#outdir<-paste("$output/",picname,sep="")
#dir.create(outdir)
#setwd(outdir)
#sample = paste(clusters,collapse=",")
#union<-union_all[,clusters]
#add end
if(length(union[,1])<=50){
	showname=TRUE;
}else{
	showname=FALSE;
}
if(length(union[1,])<=10){
	cell_widths=36
	cell_width=34
}else{
	cell_widths=floor(360/length(union[1,]))
	cell_width=floor(300/length(union[1,]))
}
num=length(union[,1])
if(dim(union)[2]==2){
	scale_row_col="column"
}else{
	scale_row_col="row"
}

if(length(union[,1])>50){
	pdf("Hcluster_heatmap.detail.pdf",height=0.015*num+8, width=6)
	pheatmap(union, color=colorRampPalette(rev(c("red","white","blue")))(100), cluster_cols=F, scale=scale_row_col,legend=T,show_rownames=TRUE, fontsize_row=1, cellwidth=cell_width, main="Cluster analysis of differentially expressed sRNA")
	dev.off()
}

pdf("Hcluster_heatmap.pdf")
pheatmap(union, color=colorRampPalette(rev(c("red","white","blue")))(100), cluster_cols=F, scale=scale_row_col,legend=T,show_rownames=showname,cellwidth=cell_widths, main="Cluster analysis of differentially expressed sRNA")
dev.off()
png("Hcluster_heatmap.png", type="cairo-png",width=480*4,height=480*4,res=72*4)
pheatmap(union, color=colorRampPalette(rev(c("red","white","blue")))(100), cluster_cols=F, scale=scale_row_col,legend=T,show_rownames=showname,cellwidth=cell_widths, main="Cluster analysis of differentially expressed sRNA")
dev.off()

# ==============================================================================
## relative expression cluster: log2(tpm ratio) cluster
# h_cluster
source("$Bin/R_hclust_log.r")
R_hclust("DE_union_for_cluster")

# SOM cluster
source("$Bin/R_som.r")
R_som("DE_union_for_cluster")

# Kmeans cluster
source("$Bin/R_kmeans.r")
R_kmeans("DE_union_for_cluster")

file.remove("Rplots.pdf")
# ==============================================================================
END

open(OUT,">$output/cluster.R");
print OUT $R;
close(OUT);

open R,"|$R_v2153/R --vanilla --slave" or die $!;
print R $R;
close R;
