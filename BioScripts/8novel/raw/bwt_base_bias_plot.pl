#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin '$Bin';
use warnings;
use Config::Tiny;

##options
my %options=();
getopts("i:a:b:W",\%options);

my $usage="usage:
\tperl $0 [options] -i read.ref.bwt
[options]
\t-a [int] max_miRNA_length during drawing first basebias map(default 30)
\t-b [int] max_miRNA_length during drawing basebias at each position(default 22)
\t-m [int] minimum length of small RNA during drawing first basebias map(default 18)
\t-W	read counts are weighted by their number of mappings. e.g. A read maps twice so each position gets 0.5 added to its read profile
\n";

if(not $options{'i'}){
	die $usage;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $R_v2153 = $Config->{srnaenv}->{R_v2153};


my $file=$options{'i'};
my $n1=30;
$n1=$options{'a'} if(defined($options{'a'}));
my $n2=22;
$n2=$options{'b'} if(defined($options{'b'}));
my $n3=18;
$n3=$options{'m'} if(defined($options{'m'}));

my %mapcounts;
if($options{'W'}){
	open(IN,"$file")||die "File $file not found";
	while(<IN>){
		if(/^(\S+)/){
			$mapcounts{$1}++;  #the same id will map ? times 
		}
	}
	close(IN);
}

my %hash=();
my %position=();
my @scores;
my $len_sc;
my ($sample, $length, $base);
open FILE,$file;  # read again
while(<FILE>){
	chomp;
# TR1_17054_x1	+	tur-mir-71-2	70	TCACTACTTTGTCTTTGGC	IIIIIIIIIIIIIIIIIII	1	
	my @line = split(/\t/);
	@scores = split(/_x/,$line[0]); #TR1_17054_
	$sample = $1 if($scores[0] =~ /^(\S+)_/);
	$len_sc = $scores[$#scores]; # 1 after x, actually "readcount"
	if($options{'W'}){$len_sc=$len_sc/$mapcounts{$line[0]};}  # above stat mapped times of the same id
	$length=length($line[4]);   # reads length
	my @bases=split("",$line[4]);  
	$base=$bases[0];
	$hash{$sample}{$length}{$base}+=$len_sc;   # first base count, %hash
	for my $i(0..$#bases){
		$position{$sample}{$i+1}{$bases[$i]}+=$len_sc;  # each base count, %position
	}

}
close FILE;

################   last read   #############################################

my $filelist1;
my $filelist2;

my @base=qw(A T C G);

###############   first base bias summary file ##############################
#############################################################################
#foreach my $sam(keys %hash)  # %hash : first
foreach my $sam(sort {$a cmp $b} keys %hash)  # %hash : first
{
	my %bias=%{$hash{$sam}};  # sam: sample, %bias fix to one sample
	open SAM,">$sam"."\.firstbase";
	$filelist1.= "$sam"."\.firstbase\t";
	print SAM "length\tA\tU\tC\tG\n";
	for my $len($n3..$n1) # $n1=30; $n3=18; $n2=22; 
	{
#	foreach my $len(sort {$a <=> $b}keys %bias){
		print SAM $len."\t";
		my @tmp;
		foreach(@base){
			if(!defined($bias{$len}{$_})){
				$bias{$len}{$_}=0;
			}
			push @tmp,$bias{$len}{$_};
		}
		my $temp=join "\t",@tmp;
		print SAM $temp."\n";
	}
	close SAM;
}

###############   base bias summary file along position##############################
#####################################################################################
my $samples;
# foreach my $sam(keys %position){
foreach my $sam(sort {$a cmp $b} keys %position)
{
	my %bias=%{ $position{$sam}};
	$samples.=$sam."\t";
	open SAM,">$sam"."\.position";
	$filelist2.="$sam"."\.position\t";
	print SAM "position\tA\tU\tC\tG\n";
	for my $pos(1..$n2) # $n2=22
	{
#	foreach my $pos(sort {$a <=> $b} keys %bias){
#		if($pos>$n){
#			last;
#		}
		print  SAM $pos."\t";
		my @tmp;
		foreach(@base){
			if(!defined($bias{$pos}{$_})) ### here, shouldn't be $pos-1 ??
			{
				$bias{$pos}{$_}=0;
			}
			push @tmp,$bias{$pos}{$_};
		}
		my $temp=join "\t",@tmp;
		print SAM $temp."\n";
	}
	close SAM;
}

my $rscript= <<"END";
###################################################################################################
files1<-unlist(strsplit("$filelist1","\\t"));
files2<-unlist(strsplit("$filelist2","\\t"));
samples<-unlist(strsplit("$samples","\\t"));

for(i in 1:length(files1)){
	firstbias<-read.table(files1[i],sep="\\t",row.names=1,header=TRUE);
	firstbias<-as.matrix(firstbias)
	firstbias_ratio<-sapply(1:nrow(firstbias),function(x) (firstbias[x,]/sum(firstbias[x,]))*100)
	colnames(firstbias_ratio)<-rownames(firstbias)
	png(paste(files1[i],"png",sep="."),type="cairo-png",height=300*5,width=(6+nrow(firstbias)*2)*25*5,res=72*5)
	barplot(firstbias_ratio,col=c("red","blue","green","darkred"),space=0.5,xlab="Length (nt)",ylab="Percent (%)",main=paste("miRNA First Nucleotide Bias ","(",samples[i],")",sep=""))
	text((1:nrow(firstbias))*1.5-0.5,100,labels=rowSums(firstbias),pos=3,xpd=TRUE,cex=1)
	legend(nrow(firstbias)*1.5,100,legend=c("G","C","U","A"),col=c("darkred","green","blue","red"),pch=15,xpd=TRUE,bty="n")	
	dev.off()
#	pdf(paste(files1[i],"pdf",sep="."))
#	barplot(firstbias_ratio,col=c("red","blue","green","darkred"),space=0.5,xlab="Length (nt)",ylab="Percent (%)",main=paste("miRNA First Nucleotide Bias ","(",samples[i],")",sep=""))
#	text((1:nrow(firstbias))*1.5-0.5,100,labels=rowSums(firstbias),pos=3,xpd=TRUE,cex=1)
#	legend(nrow(firstbias)*1.5,100,legend=c("G","C","U","A"),col=c("darkred","green","blue","red"),pch=15,xpd=TRUE,bty="n")
#	dev.off()
}


for(j in 1:length(files2)){
	position<-read.table(files2[j],sep="\\t",row.names=1,header=TRUE);
	position<-as.matrix(position)
	position_ratio<-sapply(1:nrow(position),function(x) (position[x,]/sum(position[x,]))*100)
	colnames(position_ratio)<-rownames(position)
	png(paste(files2[j],"png",sep="."),type="cairo-png",width=(6+nrow(position)*2)*15*5,height=300*5,res=72*5)
	barplot(position_ratio,col=c("red","blue","green","darkred"),space=0.5,xlab="Position",ylab="Percent (%)",main=paste("miRNA Nucleotide Bias at Each Position ","(",samples[j],")",sep=""))
#	text((1:nrow(position))*2-0.5,100,labels=rowSums(position),pos=3,xpd=TRUE,cex=0.7)
	legend(nrow(position)*1.5,100,legend=c("G","C","U","A"),col=c("darkred","green","blue","red"),pch=15,xpd=TRUE,bty="n")	
	dev.off()
#	pdf(paste(files2[j],"pdf",sep="."))
#	barplot(position_ratio,col=c("red","blue","green","darkred"),space=0.5,xlab="Position",ylab="Percent (%)",main=paste("miRNA Nucleotide Bias at Each Position ","(",samples[j],")",sep=""))
#	legend(nrow(position)*1.5,100,legend=c("G","C","U","A"),col=c("darkred","green","blue","red"),pch=15,xpd=TRUE,bty="n")
#	dev.off()
}

#####################################################################################################
END

open IN,"|$R_v2153/R --vanilla --slave" or die $!;
print IN $rscript;
close IN;

open OUT, ">base_bias.R" or die $!;
print OUT $rscript;
close OUT;
