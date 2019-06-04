#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

#plot venn, k-means cluster and heatcluster

sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description	plot venn

Options
		-in <s>	: indir
		-r <s>	: file of rowmeans of RPKM
		-compare <s>	: compare.txt
		-groupnames <s>	: group names, default group1,group2,...
		-venn <s>	: how to plot venn ; eg : inorder to get the venn of diffgeneIDs of comparision 1:2, 1:3and 2:3, and get another venn of diffgeneIDs of comparision 1:4 and 2:3, in this case need input 1:2_1:3_2:3,1:4_2:3
		-output <s>	: output dir
		-type <s>	: use which genelist eg:all/up/down. Default:all
		-h|?|help	: show this help
=========================================================================
USAGE
}
my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $R_v312 = $Config->{srnaenv}->{R_v312};

my ($in, $rowmeans, $compare, $groupnames, $venn, $output, $tp, $help);
GetOptions(
		"in=s"=>\$in,
		"r=s"=>\$rowmeans,
		"compare=s"=>\$compare,
		"groupnames=s"=>\$groupnames,
		"venn=s"=>\$venn,
		"output=s"=>\$output,
		"type=s"=>\$tp,
		"h|?|help"=>\$help,
);

if(!defined($in) || !defined($rowmeans) || !defined($compare) || !defined($groupnames) || !defined ($venn) || !defined($output) || defined($help)){
	&usage;
	exit 0;
}
$tp||= "all";
my $temp="";
if($tp eq "all"){
	$temp = "";
}elsif($tp eq "up"){
        $temp = "_up";
}elsif($tp eq "down"){
        $temp = "_down";
}

unless( -d $output){
	`mkdir -p $output`;
}

my @gn=split /,/,$groupnames;
open VC,">$output/venn.txt";
my @vc1=split /,/,$venn;
foreach(my $j=0;$j<@vc1;$j++){
	my @vc2=split /_/,$vc1[$j];
	foreach my $vc(@vc2){
		my @tmp=split /:/,$vc;
		$vc=$gn[$tmp[0]-1]."vs".$gn[$tmp[1]-1];
	}
	my $tmp1=join("_",@vc2);
	my $tmp2=join(",",@vc2);
	print VC "$tmp1\t$tmp2\n";
	&getvennlist($in,$tmp2,$output."/".$tmp1,$temp);
}
close VC;

open TYPE, ">$output/type.txt";
print TYPE "$tp\n";
close TYPE;

my $R =<< "END";
#==============================================================================
#load data
#
library("VennDiagram")
library("Vennerable")
library('reshape2')
library("ggplot2")
library("pheatmap")
rp<-read.table("$rowmeans",head=T)
#------------------------------------------------------------------------------

#==============================================================================
# prepare data for venn and cluster plotting
setwd("$in")
compare<-read.table("$compare",head=F)
array<-compare[,1]
array<-as.vector(array)
dir<-paste("$in/",array,"/",array,".Differential_analysis_results.xls",sep="")
IDdir<-paste("$in/",array,"/",array,".DElist","$temp",".txt",sep="")

un<-as.matrix("sRNA")
k<-as.matrix(rp[,1])
k<-as.matrix(k[order(k[,1]),])
colnames(k)<-"sRNA"
for(i in 1:length(array)){
	gene<-read.table(dir[i],head=T)
        n=length(gene[1,])
        gene<-gene[order(gene[,1]),]
        sig<-as.matrix(gene[,n])
	colnames(sig)<-array[i]
	k<-cbind(k,sig)
	if(is.na(file.info(IDdir[i])\$size)){
		next;
	}else{
        	a<-read.table(IDdir[i],head=F)
        	un<-rbind(un,as.matrix(a))
	}
}
un<-as.matrix(un[-1,])
index<-duplicated(un[,1])
union1<-as.matrix(un[!index,])
union<-subset(rp,sRNA %in% union1[,1])

write.table(k, file="$output/Signature.xls", row.names=F, sep="\\t", quote=F, col.names=T)

#=================================================================================
vc<-read.table("$output/venn.txt",head=F)
venn.legend <- function(venn.plot,labels,upper.ratio=0.8,fontsize=15,colors=NULL){
    if(is.null(colors))
        stop("You should specify colors")
    up <- upper.ratio
    down <- 1 - up
    pushViewport(viewport(x=0,y=down,width=1,height=up,just=c("left","bottom"),name="vennplot"))
    #grid.rect(gp=gpar(fill='red',alpha=0.2))
    grid.draw(venn.plot)
    popViewport()
    pushViewport(viewport(x=0,y=0,width=1,height=down,just=c("left",'bottom'),name='legend'))
    #grid.rect(gp=gpar(fill='blue',alpha=0.2))
    for(x in 1:length(labels)){
        grid.text(label = labels[x], x=0.1, y = unit(1,'npc') - unit(x-2.5,"lines"), just=c("left",'top'), gp=gpar(fontsize=fontsize,col=colors[x]))
    }
    popViewport()
}
for(i in 1:dim(vc)[1]){
	array<-vc[i,1]
	name<-vc[i,2]
	outdir<-paste("$output/",array,"/",array,".",sep="")
	array<-as.vector(unlist(strsplit(as.character(name),",")))
	IDdir<-paste("$in/",array,"/",array,".DElist","$temp",".txt",sep="")
	mylist=vector(mode="list",length=length(array))
	for(j in 1:length(array)){
		if(is.na(file.info(IDdir[j])\$size)){
			next
		}else{
			a<-read.table(IDdir[j],head=F)
			mylist[[j]]<-as.matrix(a)
		}
	}
#Venn diagram#####
	ven=mylist
	names(ven)<-array;
	if(length(array)==2){
        	v=venn.diagram(ven,alpha=0.4,filename=NULL,cex=1.2,cat.cex=1.2,margin=0.05,scaled=F,fill=c("Gold1", "darkorchid1"),cat.pos=c(0,0),cat.dist=c(0.02,0.02),lty=0,fontfamily="sans", cat.fontfamily="sans")
        grid.newpage()
        pdf(paste(outdir,"venn.pdf",sep=""))
        grid.draw(v)
        dev.off()
       png(paste(outdir,"venn.png",sep=""),type="cairo-png")
        grid.draw(v)
        dev.off()
	}

	if(length(array)==3){
		v=venn.diagram(ven,filename=NULL,cex=1.4,cat.cex=1.6,margin=0.05,scaled=F,fill=c("Gold1","Cyan","MediumPurple"),cat.pos=c(-15,15,180),cat.dist=rep(0.05,3),lty=0,fontfamily="sans", cat.fontfamily="sans",alpha=0.4)
        grid.newpage()
        pdf(paste(outdir,"venn.pdf",sep=""))
        grid.draw(v)
        dev.off()
        png(paste(outdir,"venn.png",sep=""),type="cairo-png")
        grid.draw(v)
        dev.off()
	}

	if(length(array)==4){
            ######## added for legend
            names(ven) <- LETTERS[1:length(array)]
            colors <- c("cornflowerblue", "green", "Gold1", "darkorchid1")
            labels <- paste(LETTERS[1:length(array)], array, sep=": ")
            ######## END 
        	 v=venn.diagram(ven, cex=1.2,cat.cex=1.4,scaled=F,filename=NULL, fontfamily="sans", cat.fontfamily="sans",  fill=colors,  margin=0.15,cat.dist=c(0.2,0.2,0.1,0.1),lty=0,alpha=0.4, cat.fontface='bold', cat.col=colors)
        grid.newpage()
        pdf(paste(outdir,"venn.pdf",sep=""))
        venn.legend(v,labels=labels,colors=colors)
        dev.off()
        png(paste(outdir,"venn.png",sep=""),type="cairo-png")
        venn.legend(v,labels=labels,colors=colors)    
        dev.off()
	}

	if(length(array)==5){
        names(ven) <- LETTERS[1:length(array)]
	colors <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
	labels <- paste(LETTERS[1:length(array)], array, sep=": ")
		v=venn.diagram(ven,cex=1.2,cat.cex=1.4,scaled=F,filename=NULL, fontfamily="sans", cat.fontfamily="sans",  fill=colors,  margin=0.15,cat.dist=c(0.2,0.25,0.2,0.2,0.25),lty=0,alpha=0.4, cat.fontface='bold', cat.col=colors)
        grid.newpage()
        pdf(paste(outdir,"venn.pdf",sep=""))
	venn.legend(v,labels=labels,colors=colors)
        dev.off()
	png(paste(outdir,"venn.png",sep=""),type="cairo-png")
	venn.legend(v,labels=labels,colors=colors)
        dev.off()
	}
}
#==============================================================================

END

open R,"|$R_v312/R --vanilla --slave" or die $!;
print R $R;
close R;

open RS,">$output/plot_venn.R";
print RS $R;
close R;

sub getvennlist{
	my($indir,$list,$out,$temp)=@_;
	unless(-d $out){
		`mkdir $out`;
	}
	my %fnhash;
	my %idhash;
	my @idfiles=split /,/,$list;
	foreach my $idlist(@idfiles){
		my $filename=$indir."/".$idlist."/".$idlist.".DElist${temp}.txt";
		open IN,"<$filename";
		while(<IN>){
			chomp;
			my $id=$_;
			push @{$idhash{$id}},$idlist;
		}
		close IN;
	}
	foreach my $ID(sort keys %idhash){
        	my $fn=join("_", @{$idhash{$ID}});
        	$fn=$fn.".venn.xls";
        	push @{$fnhash{$fn}},$ID;
	}
	foreach my $key(sort keys %fnhash){
       		open OUT,">$out/$key";
        	my $ids=join("\n",@{$fnhash{$key}});
        	print OUT "$ids\n";
        	close OUT;
	}
}
