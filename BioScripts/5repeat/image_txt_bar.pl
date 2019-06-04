#use strict;
use warnings;
use FindBin '$Bin';
use Config::Tiny;

my $input=shift;
my $o_prefix=shift;
my $ru_marker=shift;

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $R_v312 = $Config->{srnaenv}->{R_v312};

my $R= <<"END";
#========================================================================================================
library('ggplot2')
rc<-read.table("$input",sep="\\t",header=TRUE)
n<-length(colnames(rc))
labels=rc[,1]
for(i in 2:n){
	val<-round(rc[,i]/sum(rc[,i]),4)*100
	var_geom<-rc[,i]
	dat=data.frame(val,labels)
	p<-ggplot(dat,aes(x=factor(labels,levels=labels),y=val,fill=labels))+theme(panel.background=element_blank()) +theme(axis.line=element_line(colour="black"))+theme(panel.grid.major=element_line(colour=NA)) +theme(panel.grid.minor=element_line(colour=NA))+geom_bar(stat='identity',aes(width=0.8))+ylim(0,100)+labs(x='Types',y='Percent (%)',title=paste("Repeat classifiction - $ru_marker reads (",colnames(rc)[i],")"),sep="")+geom_text(aes(label=as.character(var_geom)),hjust=0,size=3)+theme(legend.position="none")+coord_flip()+scale_fill_hue(l=40)
	ggsave(filename=paste("$o_prefix","_",colnames(rc)[i],"_bar.png",sep=""),type='cairo',dpi=300)
	ggsave(filename=paste("$o_prefix","_",colnames(rc)[i],"_bar.pdf",sep=""))
dev.off()
}
#=======================================================================================================
END

open FILE, ">Rscript";
print FILE $R;
system "$R_v312/R CMD BATCH Rscript";
