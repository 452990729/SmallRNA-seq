### ARGS: [enrichment result] [output dir] [prefix].png/pdf [figure title]

### the input enrichment result like :
#GO:0009765      photosynthesis, light harvesting        biological_process      0       0       8       544     8       12000
#   comp81484_c1,comp83792_c0,comp90723_c0,comp91561_c2,comp91561_c4,comp92774_c0,comp92987_c0,comp96238_c0

args<-commandArgs(TRUE)

dat<-read.table(args[1],header=TRUE,sep="\t",fill=TRUE)

setwd(args[2])
dat<-dat[,-10]
dat<-subset(dat,dat[,5]<0.05 & dat[,6]>0)
dat<-dat[order(dat[,3]),]
ontos<-as.factor(as.character(dat[,3]))
class=levels(ontos)

sub_dat2<-function(dat,num){
        if(nrow(dat)>num){
                x<-num
        }else{
                x<-nrow(dat)
        }
        dat[1:x,]
}
dat2<-data.frame()
for(c in 1:length(class)){
        dat2<-rbind(dat2,sub_dat2(subset(dat,Term_type==class[c]),15))
}

diff<-dat2[1,7]
percent<-round(dat2[,c(6,8)]/dat2[,c(7,9)]*100,2)
out<-data.frame(GO_Term=dat2[,1],Name=dat2[,2],Ontology=dat2[,3],Percent_DEG=percent[,1],Number_DEG=dat2[,6],Percent_bg=percent[,2],Number_bg=dat2[,8])
write.table(out,file=paste(args[3],"_gene_count.txt",sep=""),quote=F,row.names=F)

title=paste("Enriched GO Terms\n","(",args[4],")",sep="");
trim_string <- function(str,num){
	  strTrim<-str
        if(nchar(str)>num){
                strTrim<-paste(substr(str,1,num),"...",sep="");
	  }
        strTrim
}
labels<-sapply(as.character(dat2[,2]),function(x) trim_string(x,35))
tick=ceiling(max(percent[,1])*1.1/10)*10
ticks=seq(0,tick,by=10)

ontos<-as.factor(as.character(dat2[,3]))
ind<-sapply(class,function(x) grep(x,as.character(ontos))[[1]])
a<-c(1,ind[-1])
b<-c(ind[-1]-1,length(ontos))
c<-rbind(a,b)
names<-class
names[names=="biological_process"]="BP"
names[names=="cellular_component"]="CC"
names[names=="molecular_function"]="MF"

colnames(c)<-names


pdf(paste(args[3],"pdf",sep="."),width=(length(percent[,1])*2.75+15)*0.2,height=25*0.2+3)

par(mar=c(40/2,10,4,5))

barplot(t(percent),space=c(0,0.75),col=c("#76BAFF","#A3FF85"),beside=TRUE,xaxt="n",yaxt="n",ylim=c(0,tick),main=title)
axis(side=1,at=(0:length(percent[,1]))*2.75+0.375,labels=rep("",length(percent[,1])+1))
axis(side=2,at=ticks,labels=ticks)
#axis(side=4,at=ticks,labels=round(ticks*diff/100))
mtext(side=2,line=2,"Percent of Genes (%)")
#mtext(side=4,line=2,"Number of Genes")
text(1:length(percent[,1])*2.75-0.25,-tick/30,labels=labels,xpd=TRUE,srt=60,pos=2)
legend(length(percent[,1]-2)*2.75,max(percent)*1.3,legend=c("DEG","Background"),col=c("#76BAFF","#A3FF85"),pch=15,border="white",bty="n",xpd=TRUE,fill=c("#76BAFF","#A3FF85"))
ln<-strwidth(labels,units="inches")*5
n<-length(labels)

o=2.75*(1:n)-1.25

q=o-ln*0.5
s=-tick/30-ln*sqrt(3)/2*tick/15
vertile<-tick/30+(max(ln)+1)*sqrt(3)/2*tick/15
t=-vertile
p=o-vertile*15/tick/sqrt(3)
tp=-tick/30-(max(ln)+1.5)*sqrt(3)/2*tick/15

for(i in 1:ncol(c)){
	start=c[1,i]
	end=c[2,i]
	segments(q[start],s[start],p[start],t,lwd=1.5,xpd=TRUE)
	segments(q[end],s[end],p[end],t,lwd=1.5,xpd=TRUE)
	segments(p[start],t,p[end],t,lwd=1.5,xpd=TRUE)
	text((p[start]+p[end])/2,tp,labels=colnames(c)[i],xpd=TRUE,pos=1)
}

dev.off()

png(paste(args[3],"png",sep="."),width=(length(percent[,1])*2.75+15)*0.2,height=25*0.2+3,units="in",type="cairo-png",res=400)

par(mar=c(40/2,10,4,5))

barplot(t(percent),space=c(0,0.75),col=c("#76BAFF","#A3FF85"),beside=TRUE,xaxt="n",yaxt="n",ylim=c(0,tick),main=title)
axis(side=1,at=(0:length(percent[,1]))*2.75+0.375,labels=rep("",length(percent[,1])+1))
axis(side=2,at=ticks,labels=ticks)
#axis(side=4,at=ticks,labels=round(ticks*diff/100))
mtext(side=2,line=2,"Percent of Genes (%)")
#mtext(side=4,line=2,"Number of Genes")
text(1:length(percent[,1])*2.75-0.25,-tick/30,labels=labels,xpd=TRUE,srt=60,pos=2)
legend(length(percent[,1]-2)*2.75,max(percent)*1.3,legend=c("DEG","Background"),col=c("#76BAFF","#A3FF85"),pch=15,border="white",bty="n",xpd=TRUE,fill=c("#76BAFF","#A3FF85"))

ln<-strwidth(labels,units="inches")*5
n<-length(labels)

o=2.75*(1:n)-1.25

q=o-ln*0.5
s=-tick/30-ln*sqrt(3)/2*tick/15
vertile<-tick/30+(max(ln)+1)*sqrt(3)/2*tick/15
t=-vertile
p=o-vertile*15/tick/sqrt(3)
tp=-tick/30-(max(ln)+1.5)*sqrt(3)/2*tick/15

for(i in 1:ncol(c)){
	start=c[1,i]
	end=c[2,i]
	segments(q[start],s[start],p[start],t,lwd=1.5,xpd=TRUE)
	segments(q[end],s[end],p[end],t,lwd=1.5,xpd=TRUE)
	segments(p[start],t,p[end],t,lwd=1.5,xpd=TRUE)
	text((p[start]+p[end])/2,tp,labels=colnames(c)[i],xpd=TRUE,pos=1)
}

dev.off()
