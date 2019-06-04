#!usr/bin/evn perl
use warnings;
use strict;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

sub usage
{
	print STDERR <<USAGE;
===============================================================================
Description     plot coExpression venn 
Version:1.0
2014-7-29      yangying\@novogene.cn

Options
	
	-venn  <s> : compares samples for venn separately,different graphs split by ',',different samples in the same graph split by ':'.(eg sample1:sample2,sample1:sample3:sample4)
	-c  <f> : cutoff of the RPKM .when RPKM more than the value then we think the gene is expressed.(default=0.3)
	-fpkm <s> :the file of fpkm (project_dir/DIFF_EXP*/Diff_analysis.out/rowmeans_fpkm.xls)
	-out <s> : output file, default ./sample_venn
===============================================================================
USAGE
}
my ($venn, $cutoff,$help,$fpkm,$out);
GetOptions(
	"h|?|help"=>\$help,
	"venn=s"=>\$venn,
	"c=f"=>\$cutoff,
        "fpkm=s"=>\$fpkm,
        "out=s"=>\$out,
);
if(!defined($venn) || !defined($fpkm) || defined($help)){
        &usage;
        exit 0;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $R_v312 = $Config->{srnaenv}->{R_v312};

my @venn=split/,/, $venn;
my $venn_num=@venn;
$cutoff ||=0.1;
my $pwd=`pwd`;
chomp($pwd);
$out ||="$pwd/sample_venn";
my $R =<< "END";
#=============================================================================d<-fpkm_filter[,1]
#                mylist[[array_name]]<-as.matrix(geneid)
#                                un<-rbind(un,as.matrix(geneid))
#                                        }
#                                        =
library("VennDiagram")
library('reshape2')
library("ggplot2")
dir.create("$out")
setwd("$out")
# prepare data for venn plotting
venn_array<-as.vector(unlist(strsplit(as.character("$venn"),",")))
for (i in 1:$venn_num){
 
	array<-as.vector(unlist(strsplit(as.character(venn_array[i]),":")))
	mylist=list()
	un<-as.matrix("geneID")
	fpkm_table<-read.delim("$fpkm")
	for(array_name in array ){
		fpkm_filter<-subset(fpkm_table,fpkm_table[[array_name]]>$cutoff)
		geneid<-fpkm_filter[,1]
		mylist[[array_name]]<-as.matrix(geneid)
		un<-rbind(un,as.matrix(geneid))
	}
	#Venn diagram#####
	ven=mylist
	venn_name<-gsub(":","_",venn_array[i])
	outdir<-paste("$out/",venn_name,sep="")
	dir.create(outdir)
	setwd(outdir)
	if(length(array)==2){
		if(length(ven[[1]]) > length(ven[[2]])) cat.pos <- c(-15,15) else cat.pos <- c(15,-15)
		v=venn.diagram(ven,alpha=0.4,filename=NULL,cex=1.2,cat.cex=1.4,margin=0.05,scaled=F,fill=c("Gold1", "darkorchid1"),cat.pos=cat.pos,cat.dist=rep(0.05,2),lty=0,fontfamily="sans", cat.fontfamily="sans")
		pdf(paste(venn_name,".DEG_Venn_diagram.pdf",sep=""))
		grid.draw(v)
		dev.off()
		png(paste(venn_name,".DEG_Venn_diagram.png",sep=""),type="cairo-png")
		grid.draw(v)
		dev.off()
		genelist=unique(c(ven[[1]],ven[[2]]))
		yes_yes<-intersect(ven[[1]],ven[[2]])
                yes_no<-ven[[1]][!(ven[[1]] %in% ven[[2]])]
                no_yes<-ven[[2]][!(ven[[2]] %in% ven[[1]])]
                write.table(as.data.frame(yes_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
                write.table(as.data.frame(genelist),paste(outdir,"/","all.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
                write.table(as.data.frame(yes_no),paste(outdir,"/",array[1],"-yes-",array[2],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
		write.table(as.data.frame(no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
}

	if(length(array)==3){
		 v=venn.diagram(ven,filename=NULL,cex=1.4,cat.cex=1.6,margin=0.05,scaled=F,fill=c("Gold1","Cyan","MediumPurple"),cat.pos=c(-15,15,180),cat.dist=rep(0.05,3),lty=0,fontfamily="sans", cat.fontfamily="sans",alpha=0.4)
		pdf(paste(venn_name,".DEG_Venn_diagram.pdf",sep=""))
		grid.draw(v)
		dev.off()
		png(paste(venn_name,".DEG_Venn_diagram.png",sep=""),type="cairo-png")
		grid.draw(v)
		dev.off()
		genelist=unique(c(ven[[1]],ven[[2]],ven[[3]]))
        	yes_yes_yes<-intersect(ven[[1]],ven[[2]])[intersect(ven[[1]],ven[[2]]) %in%  ven[[3]]]
        	yes_yes_no<-intersect(ven[[1]],ven[[2]])[!(intersect(ven[[1]],ven[[2]]) %in%  ven[[3]] )]
        	yes_no_yes<-intersect(ven[[1]],ven[[3]])[!(intersect(ven[[1]],ven[[3]]) %in%  ven[[2]])]
        	yes_no_no<-ven[[1]][!(ven[[1]] %in% c(ven[[2]],ven[[3]]))]
        	no_yes_yes<-intersect(ven[[2]],ven[[3]])[!(intersect(ven[[2]],ven[[3]]) %in% ven[[1]])]
        	no_yes_no<-ven[[2]][!(ven[[2]] %in% c(ven[[1]],ven[[3]]))]
        	no_no_yes<-ven[[3]][!(ven[[3]] %in% c(ven[[1]],ven[[2]]))]
        	write.table(as.data.frame(genelist),paste(outdir,"/","all.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_no),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_no),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_no),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	}

	if(length(array)==4){
		v=venn.diagram(ven, cex=1.2,cat.cex=1.4,scaled=F,filename=NULL, fontfamily="sans", cat.fontfamily="sans",  fill=c("cornflowerblue", "green", "Gold1", "darkorchid1"),  margin=0.15,cat.dist=c(0.2,0.2,0.1,0.1),lty=0,alpha=0.4)
		pdf(paste(venn_name,".DEG_Venn_diagram.pdf",sep=""))
		grid.draw(v)
		dev.off()
		png(paste(venn_name,".DEG_Venn_diagram.png",sep=""),type="cairo-png")
		grid.draw(v)
		dev.off()
		genelist=unique(c(ven[[1]],ven[[2]],ven[[3]],ven[[4]]))
        	yes_yes_yes_yes<-intersect(intersect(ven[[1]],ven[[2]]),intersect(ven[[3]],ven[[4]]))
        	yes_yes_yes_no<-intersect(intersect(ven[[1]],ven[[2]]),ven[[3]])[!(intersect(intersect(ven[[1]],ven[[2]]),ven[[3]]) %in% ven[[4]])]
        	yes_yes_no_yes<-intersect(intersect(ven[[1]],ven[[2]]),ven[[4]])[!(intersect(intersect(ven[[1]],ven[[2]]),ven[[4]]) %in% ven[[3]])]
        	yes_no_yes_yes<-intersect(intersect(ven[[1]],ven[[3]]),ven[[4]])[!(intersect(intersect(ven[[1]],ven[[3]]),ven[[4]]) %in% ven[[2]])]
        	no_yes_yes_yes<-intersect(intersect(ven[[2]],ven[[3]]),ven[[4]])[!(intersect(intersect(ven[[2]],ven[[3]]),ven[[4]]) %in% ven[[1]])]
        	yes_yes_no_no<-intersect(ven[[1]],ven[[2]])[!( intersect(ven[[1]],ven[[2]]) %in% c(ven[[3]],ven[[4]]))]
        	yes_no_yes_no<-intersect(ven[[1]],ven[[3]])[!( intersect(ven[[1]],ven[[3]]) %in% c(ven[[2]],ven[[4]]))]
        	no_yes_yes_no<-intersect(ven[[2]],ven[[3]])[!( intersect(ven[[2]],ven[[3]]) %in% c(ven[[1]],ven[[4]]))]
        	yes_no_no_yes<-intersect(ven[[1]],ven[[4]])[!( intersect(ven[[1]],ven[[4]]) %in% c(ven[[2]],ven[[3]]))]
        	no_yes_no_yes<-intersect(ven[[2]],ven[[4]])[!( intersect(ven[[2]],ven[[4]]) %in% c(ven[[1]],ven[[3]]))]
        	no_no_yes_yes<-intersect(ven[[3]],ven[[4]])[!( intersect(ven[[3]],ven[[4]]) %in% c(ven[[1]],ven[[2]]))]
        	yes_no_no_no<-ven[[1]][! (ven[[1]] %in% c(ven[[2]],ven[[3]],ven[[4]]))]
        	no_yes_no_no<-ven[[2]][! (ven[[2]] %in% c(ven[[1]],ven[[3]],ven[[4]]))]
        	no_no_yes_no<-ven[[3]][! (ven[[3]] %in% c(ven[[1]],ven[[2]],ven[[4]]))]
	 	no_no_no_yes<-ven[[4]][! (ven[[4]] %in% c(ven[[1]],ven[[2]],ven[[3]]))]
        	write.table(as.data.frame(genelist),paste(outdir,"/","all.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_yes_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-yes-",array[4],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_yes_no),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-yes-",array[4],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_no_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-no-",array[4],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_yes_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-yes-",array[4],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_yes_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-yes-",array[4],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_no_no),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-no-",array[4],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_yes_no),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-yes-",array[4],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_yes_no),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-yes-",array[4],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_no_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-no-",array[4],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-no-",array[4],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_no_yes_yes),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-yes-",array[4],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	 	write.table(as.data.frame(yes_no_no_no),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-no-",array[4],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_no_no),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-no-",array[4],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_no_yes_no),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-yes-",array[4],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_no_no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-no-",array[4],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")	
}

	if(length(array)==5){
		 v=venn.diagram(ven, filename=NULL, fill=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),  margin=0.15,cex=1.2,cat.cex=1.4,scaled=F,lty=0,cat.dist=c(0.2,0.25,0.2,0.2,0.25),cat.pos=c(0,-20,-160,160,20),fontfamily="sans", cat.fontfamily="sans",alpha=0.4)
		pdf(paste(venn_name,".DEG_Venn_diagram.pdf",sep=""))
		grid.draw(v)
		dev.off()
		png(paste(venn_name,".DEG_Venn_diagram.png",sep=""),type="cairo-png")
		grid.draw(v)
		dev.off()
		genelist=unique(c(ven[[1]],ven[[2]],ven[[3]],ven[[4]],ven[[5]]))
        	yes_yes_yes_yes_yes<-intersect(intersect(intersect(ven[[1]],ven[[2]]),intersect(ven[[3]],ven[[4]])),ven[[5]])
        	yes_yes_yes_yes_no<-intersect(intersect(ven[[1]],ven[[2]]),intersect(ven[[3]],ven[[4]]))[!( intersect(intersect(ven[[1]],ven[[2]]),intersect(ven[[3]],ven[[4]])) %in% ven[[5]])]
        	yes_yes_yes_no_yes<-intersect(intersect(ven[[1]],ven[[2]]),intersect(ven[[3]],ven[[5]]))[!( intersect(intersect(ven[[1]],ven[[2]]),intersect(ven[[3]],ven[[5]])) %in% ven[[4]])]
        	yes_yes_no_yes_yes<-intersect(intersect(ven[[1]],ven[[2]]),intersect(ven[[4]],ven[[5]]))[!( intersect(intersect(ven[[1]],ven[[2]]),intersect(ven[[4]],ven[[5]])) %in% ven[[3]])]
        	yes_no_yes_yes_yes<-intersect(intersect(ven[[1]],ven[[3]]),intersect(ven[[4]],ven[[5]]))[!( intersect(intersect(ven[[1]],ven[[3]]),intersect(ven[[4]],ven[[5]])) %in% ven[[2]])]
        	no_yes_yes_yes_yes<-intersect(intersect(ven[[2]],ven[[3]]),intersect(ven[[4]],ven[[5]]))[!( intersect(intersect(ven[[2]],ven[[3]]),intersect(ven[[4]],ven[[5]])) %in% ven[[1]])]
        	yes_yes_yes_no_no<-intersect(intersect(ven[[1]],ven[[2]]),ven[[3]])[! (intersect(intersect(ven[[1]],ven[[2]]),ven[[3]]) %in% c(ven[[4]],ven[[5]]) ) ]
        	yes_yes_no_yes_no<-intersect(intersect(ven[[1]],ven[[2]]),ven[[4]])[! (intersect(intersect(ven[[1]],ven[[2]]),ven[[4]]) %in% c(ven[[3]],ven[[5]]) ) ]
        	yes_no_yes_yes_no<-intersect(intersect(ven[[1]],ven[[3]]),ven[[4]])[! (intersect(intersect(ven[[1]],ven[[3]]),ven[[4]]) %in% c(ven[[2]],ven[[5]]) ) ]
        	no_yes_yes_yes_no<-intersect(intersect(ven[[2]],ven[[3]]),ven[[4]])[! (intersect(intersect(ven[[2]],ven[[3]]),ven[[4]]) %in% c(ven[[1]],ven[[5]]) ) ]
        	yes_yes_no_no_yes<-intersect(intersect(ven[[1]],ven[[2]]),ven[[5]])[! (intersect(intersect(ven[[1]],ven[[2]]),ven[[5]]) %in% c(ven[[3]],ven[[4]]) ) ]
        	yes_no_yes_no_yes<-intersect(intersect(ven[[1]],ven[[3]]),ven[[5]])[! (intersect(intersect(ven[[1]],ven[[3]]),ven[[5]]) %in% c(ven[[2]],ven[[4]]) ) ]
		no_yes_yes_no_yes<-intersect(intersect(ven[[2]],ven[[3]]),ven[[5]])[! (intersect(intersect(ven[[2]],ven[[3]]),ven[[5]]) %in% c(ven[[1]],ven[[4]]) ) ]
        	yes_no_no_yes_yes<-intersect(intersect(ven[[1]],ven[[4]]),ven[[5]])[! (intersect(intersect(ven[[1]],ven[[4]]),ven[[5]]) %in% c(ven[[2]],ven[[3]]) ) ]
        	no_yes_no_yes_yes<-intersect(intersect(ven[[2]],ven[[4]]),ven[[5]])[! (intersect(intersect(ven[[2]],ven[[4]]),ven[[5]]) %in% c(ven[[1]],ven[[3]]) ) ]
        	no_no_yes_yes_yes<-intersect(intersect(ven[[3]],ven[[4]]),ven[[5]])[! (intersect(intersect(ven[[3]],ven[[4]]),ven[[5]]) %in% c(ven[[1]],ven[[2]]) ) ]
        	yes_yes_no_no_no<-intersect(ven[[1]],ven[[2]])[! (intersect(ven[[1]],ven[[2]]) %in% c(ven[[3]],ven[[4]],ven[[5]]))]
        	yes_no_yes_no_no<-intersect(ven[[1]],ven[[3]])[! (intersect(ven[[1]],ven[[3]]) %in% c(ven[[2]],ven[[4]],ven[[5]]))]
        	no_yes_yes_no_no<-intersect(ven[[2]],ven[[3]])[! (intersect(ven[[2]],ven[[3]]) %in% c(ven[[1]],ven[[4]],ven[[5]]))]
        	yes_no_no_yes_no<-intersect(ven[[1]],ven[[4]])[! (intersect(ven[[1]],ven[[4]]) %in% c(ven[[2]],ven[[3]],ven[[5]]))]
        	no_yes_no_yes_no<-intersect(ven[[2]],ven[[4]])[! (intersect(ven[[2]],ven[[4]]) %in% c(ven[[1]],ven[[3]],ven[[5]]))]
        	no_no_yes_yes_no<-intersect(ven[[3]],ven[[4]])[! (intersect(ven[[3]],ven[[4]]) %in% c(ven[[1]],ven[[2]],ven[[5]]))]
        	yes_no_no_no_yes<-intersect(ven[[1]],ven[[5]])[! (intersect(ven[[1]],ven[[5]]) %in% c(ven[[2]],ven[[3]],ven[[4]]))]
        	no_yes_no_no_yes<-intersect(ven[[2]],ven[[5]])[! (intersect(ven[[2]],ven[[5]]) %in% c(ven[[1]],ven[[3]],ven[[4]]))]
		no_no_yes_no_yes<-intersect(ven[[3]],ven[[5]])[! (intersect(ven[[3]],ven[[5]]) %in% c(ven[[1]],ven[[2]],ven[[4]]))]
		no_no_no_yes_yes<-intersect(ven[[4]],ven[[5]])[! (intersect(ven[[4]],ven[[5]]) %in% c(ven[[1]],ven[[2]],ven[[3]]))]
        	yes_no_no_no_no<-ven[[1]][! (ven[[1]] %in% c(ven[[2]],ven[[3]],ven[[4]],ven[[5]]))]
        	no_yes_no_no_no<-ven[[2]][! (ven[[2]] %in% c(ven[[1]],ven[[3]],ven[[4]],ven[[5]]))]
        	no_no_yes_no_no<-ven[[3]][! (ven[[3]] %in% c(ven[[1]],ven[[2]],ven[[4]],ven[[5]]))]
        	no_no_no_yes_no<-ven[[4]][! (ven[[4]] %in% c(ven[[1]],ven[[2]],ven[[3]],ven[[5]]))]
	        no_no_no_no_yes<-ven[[5]][! (ven[[5]] %in% c(ven[[1]],ven[[2]],ven[[3]],ven[[4]]))]
		write.table(as.data.frame(yes_yes_yes_yes_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-yes-",array[4],"-yes-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	        write.table(as.data.frame(yes_yes_yes_yes_no),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-yes-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_yes_no_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-yes-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_no_yes_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-no-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_yes_yes_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-yes-",array[4],"-yes-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_yes_yes_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-yes-",array[4],"-yes-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_yes_no_no),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-yes-",array[4],"-no-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_no_yes_no),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-no-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	 	write.table(as.data.frame(yes_no_yes_yes_no),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-yes-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_yes_yes_no),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-yes-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_yes_no_no_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-no-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_yes_no_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-yes-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_yes_no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-yes-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_no_yes_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-no-",array[4],"-yes-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_no_yes_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-no-",array[4],"-yes-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_no_yes_yes_yes),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-yes-",array[4],"-yes-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
		write.table(as.data.frame(yes_yes_no_no_no),paste(outdir,"/",array[1],"-yes-",array[2],"-yes-",array[3],"-no-",array[4],"-no-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_yes_no_no),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-yes-",array[4],"-no-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_yes_no_no),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-yes-",array[4],"-no-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(yes_no_no_yes_no),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-no-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
       	 	write.table(as.data.frame(no_yes_no_yes_no),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-no-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_no_yes_yes_no),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-yes-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	        write.table(as.data.frame(yes_no_no_no_yes),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-no-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_no_no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-no-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
		write.table(as.data.frame(no_yes_no_no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-no-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	        write.table(as.data.frame(no_no_yes_no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-yes-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_no_no_yes_yes),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-no-",array[4],"-yes-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	        write.table(as.data.frame(yes_no_no_no_no),paste(outdir,"/",array[1],"-yes-",array[2],"-no-",array[3],"-no-",array[4],"-no-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_yes_no_no_no),paste(outdir,"/",array[1],"-no-",array[2],"-yes-",array[3],"-no-",array[4],"-no-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	        write.table(as.data.frame(no_no_yes_no_no),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-yes-",array[4],"-no-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
        	write.table(as.data.frame(no_no_no_yes_no),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-no-",array[4],"-yes-",array[5],"-no.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
	        write.table(as.data.frame(no_no_no_no_yes),paste(outdir,"/",array[1],"-no-",array[2],"-no-",array[3],"-no-",array[4],"-no-",array[5],"-yes.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

	}
}
#==============================================================================
END
open FILE,">plot_coExp_venn.Rscript";
print FILE $R;
open R,"|$R_v312/R --vanilla --slave" or die $!;
print R $R;
close R;
