R_hclust<-function(exp_data){ ###.....................
  library(NbClust)
  library(cluster)
  library(gplots)
  #library(Biobase)
  outdir = "Hclust_cluster" #creat the result file
  dir.create(outdir)
  exp_data=read.table(exp_data, header=T, sep='\t')
  #data processing and normalizing
  rownames(exp_data)=exp_data[,1]
  dims<-dim(exp_data)
  nc=dims[2]
  exp_data=exp_data[,2:nc]
  ##data normalization
  normalized_data=exp_data
  for (i in 2:dim(normalized_data)[2]){
	normalized_data[,i]=log2((exp_data[,i]+1)/(exp_data[,1]+1))
  }
  normalized_data[,1]=0
  ### find the best K value of kmeans
  #index_c=c("ball", "frey", "mcclain","kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "silhouette", "sdindex", "sdbw")
  y=NbClust(normalized_data, diss="NULL", distance = "euclidean", min.nc=4, max.nc=30, method = "ward", index = "hartigan", alphaBeale = 0.1)
  k=y$Best.nc[1,]
  if (k%%2==1){
    k=k+1
  }
  ### Generates row and column dendrograms.
  hc_sRNAs <- agnes(normalized_data, diss=FALSE, metric="euclidean", method="ward") # cluster sRNAs
  hc_samples <- hclust(as.dist(1-cor(exp_data, method="pearson")), method="complete") # cluster samples with raw data
  ### Cuts the tree and creates color vector for clusters.
  mycl <- cutree(as.hclust(hc_sRNAs), k=k)
  mycolhc <- rainbow(length(unique(mycl)), start=0.4, end=0.9)
  mycolhc <- mycolhc[as.vector(mycl)]
  names(mycolhc) <- names(mycl)
  ### creat heatmap
  nr=dim(normalized_data)[1]
  nc=dim(normalized_data)[2]
  myheatcol =rainbow(75,start=0,end=0.9)
  #feps=paste(outdir,'/h_cluster_plots.eps',sep='')
  #postscript(file=feps, horizontal=FALSE, width=8, height=18, paper="special");
  #heatmap.2(data.matrix(normalized_data), dendrogram='both', Rowv=as.dendrogram(hc_sRNAs), labRow='', Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=mycolhc, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=0.2+1/log10(nc),cexRow=0.02, margins=c(15,15), lhei=c(0.3,2), lwid=c(2.5,4))
  #dev.off()
  #fpdf=paste(outdir,'/h_cluster_plots.pdf',sep='')
  #pdf(file=fpdf, width=18, height=18, paper="special")
  #heatmap.2(data.matrix(normalized_data), dendrogram='both', Rowv=as.dendrogram(hc_sRNAs), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=mycolhc, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=0.2+1/log10(nc),labRow='', margins=c(15,15), lhei=c(0.3,2), lwid=c(2.5,4))
  #dev.off()
##########################################################################################
  ### plot the subtree and output the subclust
  plots_per_row=2
  plots_per_col=round(k/2)
  fpng=paste(outdir,'/h_cluster_plots.png',sep='')
  png(file=fpng, width=13, height=7.5*plots_per_col,units='cm', res=600, type="cairo")
  par(mfrow=c(plots_per_col, plots_per_row), pty='s')
  par(cex=0.7)
  #par(oma=c(1,1,1,1))
  #par(font.axis=2)
  for (i in 1:k){
    subcluster_i=normalized_data[mycl==i,]
    fname=paste(outdir,'/subcluster_', i, sep='')
    write.table(rownames(subcluster_i), file=fname, quote=F, sep='\t', col.names=F)
    ###plot the cluster i
    ymin=min(subcluster_i);ymax=max(subcluster_i);
    class_name= paste('subcluster_', i, sep='')
    #plot.new()
    plot_label=paste(class_name,', ', length(subcluster_i[,1]), " sRNAs", sep='')
    plot(as.numeric(subcluster_i[1,]), type='l', ylim=c(ymin, ymax), main=plot_label, col='lightgray', xaxt='n', xlab='', cex.axis=1, las=2,ylab='log2(ratio)')
    axis(side=1, at=1:length(subcluster_i[1,]), labels=colnames(subcluster_i), las=2, cex.axis=1/log10(40*nc))
    for (r in 2:length(subcluster_i[,1])){
      points(as.numeric(subcluster_i[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(subcluster_i)), type='l', col='blue')
    points(rep(0,length(subcluster_i[1,])),type='l', col='red', lwd=1)
  }
  dev.off()
  ####  pdf plot
  fpdf=paste(outdir,'/h_cluster_plots.pdf',sep='')
  pdf(file=fpdf, width=13, height=7.5*plots_per_col)
  par(mfrow=c(plots_per_col, plots_per_row), pty='s')
  par(cex=0.7)
  par(oma=c(1,1,1,1))
  for (i in 1:k){
    subcluster_i=normalized_data[mycl==i,]
    fname=paste('subcluster_', i, sep='')
    ###plot the cluster i
    ymin=min(subcluster_i);ymax=max(subcluster_i);
    plot_label=paste(fname,', ', length(subcluster_i[,1]), " sRNAs", sep='')
    plot(as.numeric(subcluster_i[1,]), type='l', ylim=c(ymin, ymax), main=plot_label, col='lightgray', cex.axis=1.2, xaxt='n', xlab='', ylab='log2(ratio)', cex.lab=1.5)
    axis(side=1, at=1:length(subcluster_i[1,]), labels=colnames(subcluster_i), las=2, cex.axis=0.7+1/log2(nc))
    for (r in 2:length(subcluster_i[,1])){
      points(as.numeric(subcluster_i[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(subcluster_i)), type='l', col='blue')
    points(rep(0,length(subcluster_i[1,])),type='l', col='red', lwd=1)
  }
  dev.off()

}
