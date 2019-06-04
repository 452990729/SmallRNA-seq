#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);
use Config::Tiny;

#my $Bin="/PUBLIC/source/RNA/smallRNA/version3/10.1diff_aly//bin/";
# Used for Plot density box, cluster 

my ($rc, $tpm, $mean_rc, $mean_tpm, $output);
GetOptions(
		"rc=s"=>\$rc,
		"tpm=s"=>\$tpm,
		"output=s"=>\$output,
);
my $usage=<<END;
------------------------------------------------------------------------------
perl $0
-rc		readcount file
-tpm		tpm file
-output		outdir
------------------------------------------------------------------------------
END
die $usage if (!$rc or !$tpm or !$output);
$rc=abs_path($rc);
$tpm=abs_path($tpm);
mkdir $output if(!-e $output);
$output=abs_path($output);
my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $R_v312 = $Config->{srnaenv}->{R_v312};

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
END

open(OUT,">$output/cluster.R");
print OUT $R;
close(OUT);

open R,"|$$R_v312  --vanilla --slave" or die $!;
print R $R;
close R;
