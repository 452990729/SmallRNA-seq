#Edit by jiangxiaoxue
#2012/10/25

use strict;
use warnings;
use FindBin '$Bin';
use Config::Tiny;

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $R_v2153= $Config->{srnaenv}->{R_v2153};
my $infile=shift;
my $outprefix=shift;

open(IN,"$infile");
my $record=0;
my @sample_mature_array;
my @family_array;
my @other_sample_array;
my %family_sample_mature;
my %family_sample_mature_count;
my %family;
my %family_count;
while(<IN>){
	chomp;
	$record++;
	my @tmp=split("\t",$_);
	if($record==1){
		@sample_mature_array=@tmp;
	}elsif($record==2){
		@family_array=@tmp;
		for my $i(1..$#tmp){
			if($tmp[$i]){
				if(!defined($family_sample_mature{$tmp[$i]})){
					$family_sample_mature_count{$tmp[$i]}=1;
					$family_sample_mature{$tmp[$i]}=$sample_mature_array[$i];
				}else{
					$family_sample_mature_count{$tmp[$i]}++;
					$family_sample_mature{$tmp[$i]}.=" ".$sample_mature_array[$i];
				}
			}
		}
	}else{
		push @other_sample_array,$tmp[0];
		for my $i(1..$#family_array){
			if($tmp[$i]){
				$family{$family_array[$i]}{$tmp[0]}=$tmp[$i];
				my @j=split(" ",$tmp[$i]);
				$family_count{$family_array[$i]}{$tmp[0]}=@j;
			}else{
				$family{$family_array[$i]}{$tmp[0]}="";
				$family_count{$family_array[$i]}{$tmp[0]}=0;
			}
		}
	}
}
close(IN);

open(OUT,">$outprefix.detail.txt");
my $family_head="";
my $family_sample_head="";
foreach my $i(keys %family_sample_mature ){
	$family_head.="\t".$i;
	$family_sample_head.="\t".$family_sample_mature{$i};
}
print OUT "$family_head\n";
print OUT "$family_sample_head\n";
foreach my $i(@other_sample_array){
	print OUT "$i";
	foreach my $j(keys %family_sample_mature ){
		print OUT "\t",$family{$j}{$i};
	}
	print OUT "\n";
}
close(OUT);


open(OUT,">$outprefix.mir_num.txt");
my $R_family_head="c(\"\"";
$family_head="";
$family_sample_head="";
foreach my $i(keys %family_sample_mature ){
	$family_head.="\t".$i;
	$R_family_head.=",\"$i\"";
	$family_sample_head.="\t".$family_sample_mature_count{$i};
}
$R_family_head.=")";
print OUT "$family_head\n";
print OUT "$family_sample_head\n";
foreach my $i(@other_sample_array){
	print OUT "$i";
	foreach my $j(keys %family_sample_mature ){
		print OUT "\t",$family_count{$j}{$i};
	}
	print OUT "\n";
}
close(OUT);

open(OUT,">$outprefix.mir_sign.txt");
$family_head="";
$family_sample_head="";
foreach my $i(keys %family_sample_mature ){
        $family_head.="\t".$i;
}
print OUT "$family_head\n";
foreach my $i(@other_sample_array){
        print OUT "$i";
        foreach my $j(keys %family_sample_mature ){
		if($family_count{$j}{$i}>0){
	                print OUT "\t+";
		}else{
			print OUT "\t-";
		}
        }
        print OUT "\n";
}
close(OUT);

my $pwd=`pwd`;
chomp($pwd);
my $R=<< "END";
#==============================================================================
library("pheatmap")
setwd("$pwd")
res<-read.table("$outprefix.mir_num.txt", sep="\\t",header=TRUE)
colnames(res)<-$R_family_head
rownames(res)<-res[,1]
res<-res[,-1]
png("$outprefix.heatmap.png", type="cairo-png",width=480*6,height=480*6,res=72*6)
pheatmap(res, cluster_cols=F,cluster_rows=T,scale="column",legend=T,show_rownames=T,show_colnames=T,main="miRNA family analysis",color=colorRampPalette(c("white","red"))(50),fontsize_row=2)
dev.off()
pdf("$outprefix.heatmap.pdf")
pheatmap(res, cluster_cols=F,cluster_rows=T,scale="column",legend=T,show_rownames=T,show_colnames=T,main="miRNA family analysis",color=colorRampPalette(c("white","red"))(50),fontsize_row=2)
dev.off()
#==============================================================================

END
;

print $R;
open R,"|$R_v2153/R --vanilla --slave" or die $!;
print R $R;
close R;
