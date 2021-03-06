#!/usr/bin/perl -w
#Perform GO enrichment analysis 
#=================================
use warnings;
use strict;
use GO::Parser;
use Cwd;
use Getopt::Long;
use File::Basename;
use Config::Tiny;

#get options from commond lines, otherwise exit
sub usage
{
  print STDERR <<USAGE;
================================================================================================================================
Description     Perform go enrichment analysis corrected the bias of gene length and plot the GO DAG graph.
				A gtf file or gene length file or sequence fasta file must be provided.
perl $0 -i <input gene_list> -goann <gene2GO annotation> -o <goseq output dir> -gtf/-length/-fasta [-enrichemntMethod "Wallenius"] [-padjust "fdr"]

         Options
		-i			Input different-gene-list file.
		-goann			gene to GO annotations file.
		-o			Output dir.
		-gtf			Gene sets for your organism, must be a gtf file.
		-length			gene length file, which the first col is gene names, and the second is length data.		
		-fasta			sequence fasta file.
		-p			prefix of the output results and figures.
		-enrichmentMethod	use to calculate category enrichment scores, "Wallenius"(default),"Sampling","Hypergeometrics".
		-padjust		P.adjust methods, "holm","hochberg","hommel","bonferroni","BH","BY","fdr(defalt)","none".
		-h|?|help		show this help.

Author: zhangmin
Revise Date: 11/28/2012
Version: goseq_v2
Changes: GO::Parser module was used to modefy the GO annotation of genes.
	 
==================================================================================================================================
USAGE
}

my ($input, $goann, $gtf, $length_file,$fasta, $help, $enrichmentM, $padjust, $output, $topcc, $topbp, $topmf, $p);

GetOptions(
	"h|?|help"	=>\$help,
	"i=s"		=>\$input,
	"o=s"		=>\$output,
	"goann=s"	=>\$goann,
	"length:s"	=>\$length_file,
	"gtf:s"		=>\$gtf,
	"fasta:s"	=>\$fasta,
	"p:s"		=> \$p,
	"enrichmentMethod:s"=>\$enrichmentM,
	"padjust:s"	=>\$padjust
);

#######################################################################################
#either $gtf or $length_file, input,$goann,$output must need, if not or help calls, exit and call help
if(defined($length_file)){
	if(!defined($input)||!defined($goann)||defined($gtf)||!defined($output)||defined($help)){
 	    &usage;
    	 exit 0;
	}
}else{
	if(!defined($gtf)||!defined($input)||!defined($goann)||!defined($output)||defined($help)){
 	    &usage;
    	 exit 0;
	}
}
if((defined($gtf) && defined($length_file)) || (defined($gtf) && defined($fasta)) || (defined($length_file) && defined($fasta))){
	print "Either gtf, length or fasta option would be provided!\n";
	&usage;
	exit 0;
}

unless (-d $output){
	!system "mkdir $output" or die "something went wrong!\n";
#	!system "cp $input $output" or warn "input file is missing!\n";
}

$enrichmentM ||="Wallenius";
$padjust ||="BH";
$p ||="ControlvsTreated";
my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $R_v312 = $Config->{srnaenv}->{R_v312};
my $GOseq = $Config->{database}->{GOseq};
my $GO = $Config->{database}->{GO};

#my $obo="$Bin/../../BioDB/GOseq/gene_ontology.1_2.obo.ab";
#my $obo_raw="$Bin/../../GO/gene_ontology.1_2.obo";
my %name=();
my %gene=();
my $gene_length;
###################################################################
######					for gtf option                   ##########
###################################################################
if(defined $gtf){
	my $bname=basename($gtf);
	$gene_length=$bname.".length";

	open(GTF,$gtf);
	open LENGTH,">$output/$gene_length";
#######calculate gene length and search for gene names
	while(<GTF>){
		chomp;
		my @array=split/\t/;
		next if ($array[2] ne "exon");
		if($array[8]=~/gene_id \"(.*?)\"/){
			if(defined($gene{$1})){
				$gene{$1}[0]=$array[3] if ($array[3]<$gene{$1}[0]);
				$gene{$1}[1]=$array[4] if ($array[4]>$gene{$1}[1]);
			}else{
				$gene{$1}[0]=$array[3];
				$gene{$1}[1]=$array[4];
			}
		}

		my ($key, $value);
       	 	if($array[8]=~/gene_id \"(.*?)\"/){
             	$key =$1;
       		}
        	if($array[8]=~/gene_name \"(.*?)\"/){
        	      $value =$1;
        	}else{
       		      $value=$key;
        	}
        	$name{$key}=$value;
	}
	close GTF;
	foreach(sort keys %gene){
		print LENGTH "$_\t";
		my $length=$gene{$_}[1]-$gene{$_}[0];
		print LENGTH "$length\n"
	}
	close LENGTH;	
}


###################################################################
######					for fasta option                 ##########
###################################################################
if(defined($fasta)){
	my $bname=basename($fasta);
	$gene_length=$bname.".length";
	open FASTA,$fasta;
	open LENGTH,">$output/$gene_length";
	
	my $header;my $sequence;
	my $firstline=1;

	while(<FASTA>){
		chomp;
		if(m/^>(.*?)\s+/){
			$header=$1;
			$name{$header}=$header;

			if(!$firstline){
				print LENGTH $header."\t".length($sequence)."\n";
				$sequence="";
			}
			$firstline=0;
		}else{
			chomp;
			$sequence.=$_;
		}
	}
	print LENGTH $header,"\t".length($sequence)."\n";
	close FASTA; close LENGTH;
}


###################################################################
######					for length option                 #########
###################################################################
if(defined($length_file)){
	my $bname=basename($length_file);
	$gene_length=$bname;
	open LN_FH,$length_file;
	open LENGTH,">$output/$gene_length";
	while(<LN_FH>){
		print LENGTH;
		chomp;
		my @array=split/\t/;
		$name{$array[0]}=$array[0];
	}
	close LN_FH; close LENGTH;
}


###############################################################################
#####################Generate a File for R Source###############################
################################################################################
open GOOUT,">$output/complete.go.txt";
my $source="Rsource";

	my $parser=new GO::Parser({handler=>'obj',use_cache=>1});
	$parser->parse($GO);
	my $graph=$parser->handler->graph;
	my $it=$graph->create_iterator;
	my %allGOs=();
	while(my $ni=$it->next_node_instance){
		my $acc=$ni->term->acc;
		my $name=$ni->term->name;
		my $onto=$ni->term->namespace();
		$allGOs{$acc}=$name."\t".$onto;
	}
	
	
	open GOANN,$goann;
	open RSOURCE,">$output/$source";

	my %go_genes=();
	my %genes_go=();
	my $goInfo="gene2GO<-list(";
	my $geneInfo="names(gene2GO)=c(";
		while(<GOANN>){
			chomp;
			my @temp=split /\t/;
			my $geneID=shift @temp;
			my %gene_GOs=();
			foreach(@temp){
				next unless defined($allGOs{$_});
				my $ref= $graph->get_reflexive_parent_terms($_);
				foreach my $term_obj(@$ref){
					my $acc=$term_obj->acc();
					$gene_GOs{$acc}=1;
					$go_genes{$acc}{$geneID}=1;
					$genes_go{$geneID}{$acc}=1
				}
			}
			my $goID = join ('","',keys %gene_GOs);
			print GOOUT $geneID."\t".join("\t",keys %gene_GOs)."\n";
			$goInfo=$goInfo."\"".$geneID."\" = c(\"".$goID."\"),";
				
			$geneInfo .="\"$geneID\",";
			
		}
	close GOANN;
	close GOOUT;
	$goInfo=~ s/,$//;
	$goInfo .=");\n";
	
	$geneInfo =~ s/,$//;
	$geneInfo .= ");\n";

	print RSOURCE "$goInfo\n$geneInfo\n";
	my $all_genes=keys %genes_go;
	close RSOURCE;



###################################################################
#######             R  => goseq; topGO                 ############
#######				go enrichemnt  =>goseq             ############
#######				go DAG graph   =>topGO			   ############
###################################################################
#
#Print the following scripts to a R file


$topcc ||= 10; $topmf ||= 10; $topbp ||= 10;

my $rOutput="$p.goseq.txt";
my $R= <<"END";
###############################################################################
#######################R 

library("goseq");
setwd("$output");

source("$source");  ##GO annotation data
universe<-names(gene2GO);
genelist <- scan("$input", what="character",quiet= TRUE);

length_data <- read.table("$gene_length",sep="\\t");
gene.len<-data.frame(len=length_data[,2])
rownames(gene.len)<-length_data[,1];
universe.len<-gene.len[universe,1]

bg_num<-length(universe);
deg_num<-sum(genelist %in% universe);

gene.vector=factor(as.integer(universe%in%genelist));
names(gene.vector)=universe;
pwf=nullp(gene.vector,bias.data=universe.len,plot.fit=FALSE);
GO=goseq(pwf,gene2cat=gene2GO,method="$enrichmentM");

#go2gene<-inverseList(gene2GO);
#bg_go<-data.frame(bg_item=sapply(1:length(go2gene),function(x) length(go2gene[[x]])))
#rownames(bg_go)<-names(go2gene);
#ind<-unlist(sapply(1:length(genelist),function(x) grep(genelist[x],names(gene2GO))))
#go2gene_sub<-inverseList(gene2GO[ind])
#deg_go<-data.frame(deg_item=sapply(1:length(go2gene_sub),function(x) length(go2gene_sub[[x]])))
#rownames(deg_go)<-names(go2gene_sub)

GO\$over_represented_pvalue[GO\$over_represented_pvalue==0]=min(GO\$over_represented_pvalue[GO\$over_represented_pvalue>0])/10000
GO\$correct<-p.adjust(GO\$over_represented_pvalue,method="$padjust");
GO\$over_represented_pvalue<- signif(GO\$over_represented_pvalue,5)
GO\$correct<-signif(GO\$correct,5)
#GO\$CAD_item<-deg_go[GO\$category,1];
#GO\$CAD_list<-rep(deg_num,nrow(GO));
#GO\$bg_item<-bg_go[GO\$category,1];
#GO\$bg_list<-rep(bg_num,nrow(GO));

write.table(GO,file="$rOutput",row.names=FALSE,quote=F,sep="\\t");
###############################################################################
END
open R,"|/PUBLIC/software/RNA/R-3.1.2/R-3.1.2/bin/R --vanilla --slave" or die $!;
print R $R;
close R;
open R,">$output/goseq.R";
print R $R;
close R;

my %go_diffgenes=();
my $diff_genes=0;
open LIST,$input;
while(<LIST>){
	chomp;
	my $id=$_;
	if(exists($genes_go{$id})){
		$diff_genes++;
		foreach(keys %{$genes_go{$id}}){
			$go_diffgenes{$_}{$id}=1;
		}
	}	
}
close LIST;

open GOSEQ,"$output/$rOutput";
<GOSEQ>;
my $out="$p.DEG_GO_enrichment_result.xls";
open OUTPUT, ">$output/$out";
print OUTPUT "GO_accession\tDescription\tTerm_type\tOver_represented_pValue\tCorrected_pValue\tDEG_item\tDEG_list\tBg_item\tBg_list\tGene_names\n";
while(<GOSEQ>){
	chomp;
#	my ($acc,$over,$under,$padjust) = split /\s+/;edit by jc in 2015/5/19
	my ($acc,$over,$under,$padjust) = (split /\t/)[0,1,2,7];
	my $term = $allGOs{$acc};
	next unless($term);
	my $diff_terms = keys %{ $go_diffgenes{$acc}};
	my $all_terms = keys %{ $go_genes{$acc}};
	if($diff_terms == 0){
		next;
	}
	print OUTPUT $acc."\t".$term."\t".$over."\t".$padjust."\t".$diff_terms."\t".$diff_genes."\t".$all_terms."\t".$all_genes."\t".join(",",keys %{ $go_diffgenes{$acc}})."\n";
}
close GOSEQ;
close OUTPUT;


my $R_topGO= <<"END";
##===================================================================================
setwd("$output")
library(GO.db)
library(goTools)
library(topGO)
library(goseq)

goseq<-read.table("$out",header=TRUE,sep="\\t")

source("$source");  ##GO annotation data
universe<-names(gene2GO);
genelist <- scan("$input", what="character",quiet= TRUE);
gene.vector=factor(as.integer(universe %in% genelist));
names(gene.vector)=universe;

goseq<-subset(goseq,Corrected_pValue<0.05)
rownames(goseq)<-goseq\$GO_accession;

if(sum(goseq\$Term_type=="biological_process")>0){
	GObpdata<-new("topGOdata",description="BP",ontology="BP",allGenes=gene.vector,annot=annFUN.gene2GO,gene2GO=gene2GO,nodeSize=10);
	bpNodes<-nodes(graph(GObpdata));
	bpScores<-goseq[bpNodes,5];
	names(bpScores)<-bpNodes;
	if(sum(!is.na(bpScores))>$topbp){
		bps<-$topbp;
	}else{
		bps<-sum(!is.na(bpScores));
	}
	if(bps>=1){
		pdf("$p.DEG_Enriched_GO_bp_DAG.pdf");
		showSigOfNodes(GObpdata,useInfo="all",bpScores,firstSigNodes=bps);
		dev.off();

		png("$p.DEG_Enriched_GO_bp_DAG.png",type="cairo-png",width=480*4,height=480*4,res=72*4);
		showSigOfNodes(GObpdata,useInfo="all",bpScores,firstSigNodes=bps);
		dev.off();
	}
}

if(sum(goseq\$Term_type=="cellular_component")>0){
	GOccdata<-new("topGOdata",description="CC",ontology="CC",allGenes=gene.vector,annot=annFUN.gene2GO,gene2GO=gene2GO,nodeSize=10);
	ccNodes<-nodes(graph(GOccdata));
	ccScores<-goseq[ccNodes,5];
	names(ccScores)<-ccNodes;
	if(sum(!is.na(ccScores))>$topcc){
		ccs<-$topcc;
	}else{
		ccs<-sum(!is.na(ccScores));
	}
	if(ccs>=1){
		pdf("$p.DEG_Enriched_GO_cc_DAG.pdf");
		showSigOfNodes(GOccdata,useInfo="all",ccScores,firstSigNodes=ccs);
		dev.off();

		png("$p.DEG_Enriched_GO_cc_DAG.png",type="cairo-png",width=480*4,height=480*4,res=72*4);
		showSigOfNodes(GOccdata,useInfo="all",ccScores,firstSigNodes=ccs);
		dev.off();
	}
}
if(sum(goseq\$Term_type=="molecular_function")>0){
	GOmfdata<-new("topGOdata",description="MF",ontology="MF",allGenes=gene.vector,annot=annFUN.gene2GO,gene2GO=gene2GO,nodeSize=10);
	mfNodes<-nodes(graph(GOmfdata));
	mfScores<-goseq[mfNodes,5];
	names(mfScores)<-mfNodes;
	if(sum(!is.na(mfScores))>$topmf){
		mfs<-$topmf;
	}else{
		mfs<-sum(!is.na(mfScores));
	}
	if(mfs>=1){
		pdf("$p.DEG_Enriched_GO_mf_DAG.pdf");
		showSigOfNodes(GOmfdata,useInfo="all",mfScores,firstSigNodes=mfs);
		dev.off();
	
		png("$p.DEG_Enriched_GO_mf_DAG.png",type="cairo-png",width=480*4,height=480*4,res=72*4);
		showSigOfNodes(GOmfdata,useInfo="all",mfScores,firstSigNodes=mfs);
		dev.off();
	}
}

#======================================================================================
END
open RBAR,"|$R_v312/R --vanilla --slave" or die $!;
print RBAR $R_topGO;
close RBAR;

open R,">$output/topGO.R";
print R $R_topGO;
close R;

my $rscript="$R_v312/Rscript";
my $bar1_r="$Bin/goseq_bar_plot.R";
my $bar2_r="$Bin/goseq_bar_plot2.R";

open R, ">$output/Enriched_GO_bar.R";
print R "#----------------- Real run use this one command ----------------------------------\n";
print R "$rscript $bar1_r $output/$out $output $p.DEG_Enriched_GO_classification $p\n\n";
print R "#----------------- Following R only for check and study, NOT run ----------------------------------\n";
print R `more $bar1_r`;
close R;

system("$rscript $bar1_r $output/$out $output $p.DEG_Enriched_GO_classification $p");
