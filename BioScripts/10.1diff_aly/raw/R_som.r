################# usage of this function   ###################
### parameter "exp_data" can be the absolute path of the data including the data name or only be the data name in the current directory
### the data formation should be like below:
###   sRNAid sample1 sample2 ...
###   sRNA1  124     321 ...
###   sRNA2  222     342 ...
### when you start R, run the function like this:
### > source('R_som.r')
### > R_som('combine.tpm.txt')
R_som<-function(exp_data){     
  library(som)
  library(gplots)
  outdir=  "SOM_cluster"#creat the result file
  dir.create(outdir)
  exp_data=read.table(exp_data, header=T, sep='\t')
  #data processing and  normalizing
  rownames(exp_data)=exp_data[,1]
  dims<-dim(exp_data)
  nc=dims[2]
  exp_data=exp_data[,2:nc]
  #if (nc>3){
   # normalized_data=normalize(exp_data, byrow=TRUE)
    #normalized_data = t(scale(t(exp_data), center=T, scale=T)) ##not suit to data with 2 samples      
  #}
  #else 
  normalized_data=exp_data
  for (i in 2:dim(normalized_data)[2]){
	normalized_data[,i]=log2((normalized_data[,i]+0.00001)/(normalized_data[,1]+0.00001))
  }
  normalized_data[,1]=0		 
  # som training
  data_som=som(normalized_data,xdim=5, ydim=6, init="linear", alpha=NULL, alphaType="inverse", neigh="gaussian", topol="rect")
  som_num=data_som$code.sum
  class_num=som_num[som_num[,3]!=0,]
  index=data_som$visual
  index_x=index[,1]
  index_y=index[,2]
  num=dim(class_num)[1] ###number of cluster
  # fetching clusters and plotting them
  plots_per_row=2
  plots_per_col=round(num/2)
  fpng=paste(outdir,'/SOM_cluster.png',sep='')
  png(file=fpng, width=13, height=6.5*plots_per_col,units='cm', res=600, type="cairo")
  par(mfrow=c(plots_per_col, plots_per_row), pty='m')
  #par(cex.axis=0.6, cex.lab=0.7, cex.sub=0.8)
  for (i in 1:num){
    x=index_x==class_num[i,1]
    y=index_y==class_num[i,2]
    subcluster_i=normalized_data[x*y==1,]
    fname=paste(outdir,'/subcluster_', class_num[i,1]+1,'_',class_num[i,2]+1, sep='')
    write.table(subcluster_i, file=fname, quote=F, sep='\t', col.names=F)
  ###plot the cluster i
    ymin=min(subcluster_i);ymax=max(subcluster_i);
    class_name= paste('subcluster_', class_num[i,1]+1,'_',class_num[i,2]+1, sep='')
    plot_label=paste(class_name,', ', length(subcluster_i[,1]), " sRNAs", sep='')
    plot(as.numeric(subcluster_i[1,]), type='l', cex.axis=1, ylim=c(ymin, ymax), las=2, main=plot_label, col='lightgray', xaxt='n', xlab='', ylab='log2(ratio)')
    axis(side=1, at=1:length(subcluster_i[1,]), labels=colnames(subcluster_i),cex.axis=1)
    for (r in 2:length(subcluster_i[,1])){
      points(as.numeric(subcluster_i[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(subcluster_i)), type='o', col='blue')
	points(rep(0,length(subcluster_i[1,])),type='l', col='red', lwd=1)
}
  dev.off()
  
  ####  pdf plot
  fpdf=paste(outdir,'/SOM_cluster.pdf',sep='') 
  pdf(file=fpdf, width=13, height=6.5*plots_per_col)
  par(mfrow=c(plots_per_col, plots_per_row), pty='s')
  #par(cex=0.7)
  #par(cex.axis=0.6, cex.lab=0.7, cex.sub=0.8)
  #par(oma=c(1,1,1,1))
  for (i in 1:num){
    x=index_x==class_num[i,1]
    y=index_y==class_num[i,2]
    subcluster_i=normalized_data[x*y==1,]
    fname=paste('subcluster_', class_num[i,1]+1,'_',class_num[i,2]+1, sep='')
    ###plot the cluster i
    ymin=min(subcluster_i);ymax=max(subcluster_i);
    plot_label=paste(fname,', ', length(subcluster_i[,1]), " sRNAs", sep='')
    plot(as.numeric(subcluster_i[1,]), type='l', cex.axis=1.7,ylim=c(ymin, ymax), las=2, main=plot_label, cex.main=2.5, col='lightgray', xaxt='n', xlab='', ylab='log2(ratio)', cex.lab=2)
    axis(side=1, at=1:length(subcluster_i[1,]), labels=colnames(subcluster_i), cex.axis=1.7)
    for (r in 2:length(subcluster_i[,1])){
      points(as.numeric(subcluster_i[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(subcluster_i)), type='o', col='blue')
	points(rep(0,length(subcluster_i[1,])),type='l', col='red', lwd=1)
  }
  dev.off()
}
##########################################################
#############    author: lili    #########################
#############    time: 2012.0926   #########################
