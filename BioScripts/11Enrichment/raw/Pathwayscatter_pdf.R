args <- commandArgs(TRUE);##sampleName ##output dir ##identify file
pdf(paste(args[2],"/",args[1],".","DEG_enriched_KEGG_pathway_scatterplot.pdf",sep=""),height=6,width=8)
library(ggplot2)
path=read.table(args[3],sep="\t",header=T)
rownames(path)<-path[,1]
colnames(path)<-c("Pathway_term","Rich_factor","Qvalue","Gene_number")
p <- ggplot(path, aes(Rich_factor,Pathway_term))
p + geom_point(aes(colour=Qvalue,size=Gene_number))+scale_colour_gradientn(colours=rainbow(3),guide = "colourbar") + labs(title="Statistics of Pathway Enrichment") + xlab("Rich factor") + ylab("Pathway term")
dev.off()
