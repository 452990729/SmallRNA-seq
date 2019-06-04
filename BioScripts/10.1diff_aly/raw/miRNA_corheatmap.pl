#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use Config::Tiny;

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $R_v2153 = $Config->{srnaenv}->{R_v2153};
my ($help, $outdir, $association, $type, $corheat, $novenn, $family, $ppi, $mirna, $mrna, $pairs, $ann, $ko, $go, $len, $genome, $kegg, $gtf, $compare, $tpm,  $fullname, $list, $enrichdir, $geneid, $genefa);
GetOptions(
	"h|help"	=>\$help,
	"tpm=s"		=>\$tpm,
);
my $usage=<<END;
-----------------------------------------------------------------------------------------
perl $0
-h|help		help
-tpm=s		miRNA TPM file
-----------------------------------------------------------------------------------------
END
die $usage if ($help or !$tpm);
my $wkdir=getcwd();
open Corheat, ">$wkdir/miRNA_TPM_corheat.R";
my $corheat_R=<<END;
# ---------------------------------------------------------------------------------------
setwd("$wkdir")
library("ggplot2");
library("reshape2");
tpm=read.table("$tpm", header=T, row.names=1);
#cor_table=cor(tpm);  cor_table
cor_table=cor(log10(tpm+1))^2;  cor_table
dat=melt(cor_table);
colnames(dat)=c("id1", "id2", "corr");
sam_num=ncol(tpm);

if(sam_num<5){
    size_number=5
}else if(sam_num>=5 && sam_num<10){
    size_number=4
}else if(sam_num>=10 && sam_num<15){
    size_number=3
}else if(sam_num>=15 && sam_num<18){
    size_number=2
}else{
    size_number=1.5
}

p=ggplot(dat, aes(id1, id2)) + geom_tile(aes(fill=corr), colour="black") + scale_fill_gradient(name=expression(R^2), low="white", high="#4876FF") + labs(title="Peason correlation between samples", x="", y="")  + theme(axis.text.x=element_text(angle=45)) + coord_fixed() +geom_text(aes(label=round(corr,3)), size=size_number);
p;
ggsave("miRNA_TPM_corheat.png", type="cairo")
ggsave("miRNA_TPM_corheat.pdf")
file.remove("Rplots.pdf")
# ---------------------------------------------------------------------------------------
END
print Corheat "$corheat_R";

chdir "$wkdir";
`$R_v2153/Rscript miRNA_TPM_corheat.R`;

