################# usage of this function   ###################
### parameter "exp_data" can be the absolute path of the data including the data name or only be the data name in the current directory
### the data formation should be like below:
###   sRNAid sample1 sample2 ...
###   sRNA1  124     321 ...
###   sRNA2  222	   342 ...
### when you start R, run the function like this:
### > source('R_kmeans.r')
### > R_kmeans('combine.tpm.txt')
R_kmeans<-function(exp_data){     ###当前文件夹下的数据名称，或此数据的绝对路径
  library(NbClust)
  outdir=  "K_means_cluster"#creat the result file
  exp_data=read.table(exp_data, header=T, sep='\t')
  #data processing and  correlation
  rownames(exp_data)=exp_data[,1]
  dims<-dim(exp_data)
  nc=dims[2]
  exp_data=exp_data[,2:nc]
  #if (nc>3){
   # normalized_data = t(scale(t(exp_data), center=T, scale=T)) ##not suit to data with 2 samples    
  #}
  #else normalized_data=exp_data
  normalized_data=exp_data
  for (i in 2:dim(normalized_data)[2]){
	normalized_data[,i]=log2((normalized_data[,i]+0.00001)/(normalized_data[,1]+0.00001))
  }
  normalized_data[,1]=0			
  ### find the best K value of kmeans
  #index_all=c("mcclain","kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "rubin")
#  y=NbClust(normalized_data, diss=diss, distance ="NULL", min.nc=2, max.nc=30, method = "kmeans", index =index_all, alphaBeale = 0.1)
  if (dims[1]>25){
    maxn=25
  }
  else maxn=round(dims[1]/2)
  y=NbClust(normalized_data, diss="NULL", distance ="euclidean", min.nc=2, max.nc=maxn, method = "kmeans", index ="hartigan", alphaBeale = 0.1)
  k=y$Best.nc[1,]
  if (k%%2==1)
	k=k+1		
  cl=kmeans(normalized_data, k, iter.max=50)
  index=cl$cluster
  plots_per_row=2
  plots_per_col=round(k/2)
  dir.create(outdir)
  fpng=paste(outdir,'/K_means_cluster.png',sep='')
  png(file=fpng, width=13, height=6.5*plots_per_col,units='cm', res=600, type='cairo')
  par(mfrow=c(plots_per_col, plots_per_row), pty='m')
  #par(cex=0.7)
  #par(oma=c(1,1,1,1))
  #par(font.axis=2)
  for (i in 1:k){
    subcluster_i=normalized_data[index==i,]
    fname=paste(outdir,'/subcluster_', i, sep='')
    write.table(subcluster_i, file=fname, quote=F, sep='\t', col.names=F)
    ###plot the cluster i
    ymin=min(subcluster_i);ymax=max(subcluster_i);
    class_name= paste('subcluster_', i, sep='')
    plot_label=paste(class_name,', ', length(subcluster_i[,1]), " sRNAs", sep='')
    plot(as.numeric(subcluster_i[1,]), type='l', ylim=c(ymin, ymax),cex.axis=1, las=2, main=plot_label, col='lightgray', xaxt='n', xlab='', ylab='log2(ratio)')
    axis(side=1, at=1:length(subcluster_i[1,]), labels=colnames(subcluster_i), cex.axis=1)
    #x=c(1,2)
    #text(x, -10, srt=60, xpd=T, pos=2)
    for (r in 2:length(subcluster_i[,1])){
      points(as.numeric(subcluster_i[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(subcluster_i)), type='o', col='blue')
	points(rep(0,length(subcluster_i[1,])),type='l', col='red', lwd=1)
  }
  dev.off()
  ####  pdf plot
  fpdf=paste(outdir,'/K_means_cluster.pdf',sep='') 
  pdf(file=fpdf, width=13, height=6.5*plots_per_col)
  par(mfrow=c(plots_per_col, plots_per_row), pty='s')
  #par(cex=0.7)
  #par(oma=c(1,1,1,1))
  for (i in 1:k){
    subcluster_i=normalized_data[index==i,]
    fname=paste('subcluster_', i, sep='')
    ###plot the cluster i
    ymin=min(subcluster_i);ymax=max(subcluster_i);
    plot_label=paste(fname,', ', length(subcluster_i[,1]), " sRNAs", sep='')
    plot(as.numeric(subcluster_i[1,]), type='l', ylim=c(ymin, ymax), cex.axis=1.7, las=2, main=plot_label, cex.main=2.5, col='lightgray', xaxt='n', xlab='', ylab='log2(ratio)', cex.lab=2)
    axis(side=1, at=1:length(subcluster_i[1,]), labels=colnames(subcluster_i), cex.axis=1.7)
    for (r in 2:length(subcluster_i[,1])){
      points(as.numeric(subcluster_i[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(subcluster_i)), type='o', col='blue')
	points(rep(0,length(subcluster_i[1,])),type='l', col='red', lwd=1)

  }
  dev.off()

##########################################################
#############    author: lili    #########################
#############    time: 2012.08   #########################
}
