#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use Config::Tiny;
use FindBin '$Bin';

my ($help, $fa, $outdir, $prefix, $min, $max);
GetOptions(
	"h|help" =>\$help,
	"fa=s"	=>\$fa,
	"od=s"	=>\$outdir,
	"pre=s"	=>\$prefix,
	"min=i"	=>\$min,
	"max=i"	=>\$max,
);
my $usage=<<END;
--------------------------------------------------------
perl $0
-h|help		help
-fa=s		fasta file
-od=s		outdir
-pre=s		prefix of outfile, eg: sample name
-min=i		min length need 
-max=i		max length need 
--------------------------------------------------------
END
die $usage if ($help or !$fa or !$outdir or !$prefix or !$min or !$max);
mkdir $outdir if (!-e $outdir);
$outdir=abs_path($outdir);

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $R_v2153 = $Config->{srnaenv}->{R_v2153};


my $pwd = `pwd`;
$pwd =~ s/\s+//g;
open F, "$outdir/$fa" or die $!;
open O, ">$outdir/$prefix\_len.stat" or die $!;
print O "Length\tCount\n";
my %len;
while(<F>)
{
	chomp;
	if(/>(\S+)\((\d+)/)
	{
		my $seqlen=length $1;
		$len{$seqlen}+=$2;
	}
}
close F;

for ($min..$max)
{
	if(!$len{$_})
	{
		$len{$_}=0;
	}
	print O "$_\t$len{$_}\n";
}
close F;
close O;

my $R=<<END;
# ------------------------------------------------------------------
setwd("$outdir")
library("ggplot2")
library("scales");
dat=read.table("$prefix\_len.stat", header=T);
Ratio=round(100*dat\$Count/sum(dat\$Count), 2)
dat\$Ratio=Ratio;

g=ggplot(dat, aes(x=factor(Length), y=Count/sum(Count), fill=factor(Length))) + geom_bar(stat="identity") +theme(legend.position="none")+theme(panel.background=element_blank()) +theme(axis.line=element_line(colour="black")) + scale_y_continuous(labels=percent) + geom_text(aes(label=Ratio), size=2, vjust=-0.5) + labs(x="Length(nt)", y="Frequence percent", title="Sequence length distribution\\n($prefix)")
ggsave(plot=g, file="$outdir/$prefix\_seq_len_distribution.png", type="cairo")
ggsave(plot=g, file="$outdir/$prefix\_seq_len_distribution.pdf")
file.remove("$outdir/Rplots.pdf")
# ------------------------------------------------------------------
END

open R, ">$outdir/len_dis_plot.R" or die $!;
print R "$R";
close R;

`$R_v2153/Rscript $outdir/len_dis_plot.R`;
