###plot the correlation between two samples
 corr_plot<-function(exp_data){
  library(ggplot2)
 library(reshape2)
  dirname<-unlist(strsplit(exp_data,basename(exp_data)))
  outdir= paste(dirname, "corr_plot",sep='')#creat the result file
  dir.create(outdir)
  exp_data=read.table(exp_data, header=T, sep='\t')
  final_name = names(exp_data)[-1] # add
  rownames(exp_data)=exp_data[,1]
  dims<-dim(exp_data)
  nc=dims[2]
   sam_num=nc-1
  #exp_data=data.matrix(exp_data[,2:nc])
  #exp_data=exp_data[rowSums(exp_data)!=0,]
  exp_data=exp_data[,2:nc]
  dat_all=log10(exp_data+1)
  sam_num=nc-1
  plot_num=sam_num*(sam_num-1)/2
  cor_table1=matrix(1,sam_num,sam_num)
  cor_table2=matrix(1,sam_num,sam_num)
  cor_table3=matrix(1,sam_num,sam_num)
  #colnames(cor_table)=rownames(cor_table)=colnames(dat_all)
  for (i in 1:(sam_num-1)){
    for (j in (i+1):sam_num){
      dat=dat_all[,c(i,j)]
      loc=max(dat[,2]) 
      dat=data.frame(dat)
      p<- ggplot(dat)
      p<- p + aes_string(x=colnames(dat)[1],y=colnames(dat)[2])
      model <- coef(lm(dat[[2]] ~ dat[[1]], data = dat))
      intercept_val <- as.numeric(model)[1]
      slope_val <- as.numeric(model)[2]
      p<- p + geom_point(size=1.5,alpha=0.3,colour="#4876FF") + geom_abline(intercept=intercept_val,slope=slope_val,linetype=2,colour="#FF7F50") + geom_rug(size=0.5,alpha=0.01,colour="#4876FF")
      R2 <- signif((cor(dat[ ,1],dat[ ,2],method="pearson"))^2,3)
      cor_table1[i,j]=R2
      cor_table1[j,i]=R2
      rho <- signif(cor(dat[ ,1],dat[ ,2],method="spearman"),3)
      cor_table2[i,j]=rho
      cor_table2[j,i]=rho
      tao <- signif(cor(dat[ ,1],dat[ ,2],method="kendall"),3)
      cor_table3[i,j]=tao
      cor_table3[j,i]=tao
      p <- p + labs(title=paste(colnames(dat_all)[i]," vs ",colnames(dat_all)[j]))
      p <- p + xlab(paste("log10(TPM+1),"," (",colnames(dat_all)[i],")")) + ylab(paste("log10(TPM+1),"," (",colnames(dat_all)[j],")"))
      y1=loc
      y2=loc-0.2
      y3=loc-0.4
      p<- p+ annotate("text",adj=0,x=0.02,y=y1, label=paste("R^2==",R2), parse=TRUE)
      p<- p+ annotate("text",adj=0,x=0.02,y=y3, label=paste("Tau=",tao))
      p<- p+ annotate("text",adj=0,x=0.02,y=y2, label=paste("Rho=",rho))
	
	p<- p+ theme(panel.background=element_blank())
	p<- p+ theme(axis.line=element_line(colour="black"))
	p<- p+ theme(legend.position=c(0.9,0.85))
	p<- p+ theme(panel.grid.major=element_line(colour=NA))
	p<- p+ theme(panel.grid.minor=element_line(colour=NA))	
      fpdf=paste(outdir,'/',colnames(dat_all)[i],'_vs_',colnames(dat_all)[j],'.scatter.pdf',sep='')
      fpng=paste(outdir,'/',colnames(dat_all)[i],'_vs_',colnames(dat_all)[j],'.scatter.png',sep='')
      ggsave(filename=fpdf, plot=p)
      ggsave(filename=fpng,type="cairo-png", plot=p)
    }
  }
 cor_table1<-data.frame(cor_table1)
 colnames(cor_table1)<-colnames(dat_all)
 cor_table1$coefficient<-colnames(dat_all)
 cor_table1<-cor_table1[,c(sam_num+1,1:sam_num)]
 names(cor_table1)[1]<-"R^2"
 ft=paste(outdir,'/cor_pearson.xls',sep='')
 write.table(cor_table1,file=ft,quote=F,row.name=F, sep="\t")
 cor_table2<-data.frame(cor_table2)
 colnames(cor_table2)<-colnames(dat_all)
 cor_table2$coefficient<-colnames(dat_all)
 cor_table2<-cor_table2[,c(sam_num+1,1:sam_num)]
 names(cor_table2)[1]<-"Rho"
 ft2=paste(outdir,'/cor_spearman.xls',sep='')
 write.table(cor_table2,file=ft2,quote=F,row.name=F, sep="\t")
 cor_table3<-data.frame(cor_table3)
 colnames(cor_table3)<-colnames(dat_all)
 cor_table3$coefficient<-colnames(dat_all)
 cor_table3<-cor_table3[,c(sam_num+1,1:sam_num)]
 names(cor_table3)[1]<-"Tau"
 ft3=paste(outdir,'/cor_kendall.xls',sep='')
 write.table(cor_table3,file=ft3,quote=F,row.name=F, sep="\t")
 if(sam_num<=6){
        size_number=5
}else if(sam_num<10){
        size_number=3
}else if(sam_num<15){
        size_number=2.5
}else if(sam_num<18){
        size_number=2
}else{
        size_number=1.5
}
heat<-cor_table1
order<-heat[,1]
order<-as.vector(as.character(order))
df<-melt(heat)
colnames(df)<-c("sample1","sample2","correlation")
df$sample1 <- factor(df$sample1,levels=final_name) # add
p<-ggplot(df,aes(sample1,sample2,label=correlation))+
geom_tile(aes(fill = correlation),colour="white") +
scale_fill_gradient(name=expression(R^2),low="white",high="#4876FF")+
theme(panel.background = element_rect(fill='white', colour='white')) +
labs(x="",y="", title="Pearson correlation between samples")+
theme(legend.position="right",axis.text.x=element_text(angle=45,vjust=1,hjust=1))+coord_fixed()+
geom_text(size=size_number)+xlim(order)+ylim(order)
ggsave(filename=paste(outdir,'/cor_pearson.pdf',sep=''),plot=p, height=max(6,round(nrow(sam_num)/3)), width=max(6,round(nrow(sam_num)/3)))
ggsave(filename=paste(outdir,'/cor_pearson.png',sep=''),type="cairo-png", plot=p, height=max(6,round(nrow(sam_num)/3)), width=max(6,round(nrow(sam_num)/3)))

}

## debug(corr_plot)
## corr_plot("/TJPROJ1/RNA/zhangjie/NHR140179_dashu_refsRNA/sRNA/11.diff/rno/merged.tpm")
