#!/usr/bin/perl -w
use warnings;
use strict;
use Cwd;
use File::Basename;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

my ($diff, $ko, $help, $method, $padjust, $out, $cutoff, $group_name, $ptype);

GetOptions(
	"diff=s"	=> \$diff,
	"ko=s"		=> \$ko,
	"help|h|?"	=> \$help,
	"m=s"		=> \$method,
	"c=i"		=> \$cutoff,
	"padj=s"	=> \$padjust,
	"out=s"		=> \$out,
	"g=s"		=> \$group_name,
	"t=s"		=> \$ptype
);

if(!defined($diff) || !defined($ko) || defined($help)){
	&usage;
	exit 0;
}

$out ||=getcwd();
$method ||="h";
$padjust ||="BH";
$cutoff ||="5";
$group_name ||=basename($diff);

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $python_v2710 = $Config->{srnaenv}->{python_v2710};
my $R_v312 = $Config->{srnaenv}->{R_v312};
my $Rscript="$R_v312/Rscript";

my $ko_extract="$Bin/ko_extract.pl";
my $identify="$Bin/noref_identify.py";
my $add_pl="$Bin/noref_add_id.pl";
my $top20_pl="$Bin/noref_extract_top20.pl";
my $plot="$Bin/Pathwayscatter.R";


unless(-d $out){
	!system("mkdir $out") or warn "Something goes wrong:$!!\n";
}

my $fg="$group_name.fg";
my $kegg_enrich="$group_name.DEG_KEGG_pathway_enrichment_result.xls";
my $kegg_enrich_add="$group_name.DEG_KEGG_pathway_enrichment_add.xls";
my $enrich_top20="$group_name.DEG_enriched_KEGG_pathway_top20.xls";


print "$perlExec ${ko_extract} $diff $ko $out/$fg\n";
print "${python_v2710} ${identify} -f $out/$fg  -b $ko -d K -m $method -n $padjust -o $out/$kegg_enrich -c $cutoff -t $ptype\n";
print "$perlExec $add_pl $ko $out/$kegg_enrich > $out/$kegg_enrich_add\n";
print "$perlExec $top20_pl $out/$kegg_enrich_add > $out/$enrich_top20\n";
print "${Rscript} $plot $group_name $out $out/$enrich_top20\n";

!system("$perlExec ${ko_extract} $diff $ko $out/$fg") or die $!;
!system("${python_v2710} ${identify} -f $out/$fg  -b $ko -d K -m $method -n $padjust -o $out/$kegg_enrich  -c $cutoff -t $ptype") or die $!;
!system("$perlExec $add_pl $ko $out/$kegg_enrich > $out/$kegg_enrich_add") or die $!;
!system("$perlExec $top20_pl $out/$kegg_enrich_add > $out/$enrich_top20") or die $!;

!system("${Rscript} $plot $group_name $out $out/$enrich_top20") or die $!;


sub usage{
	print STDERR <<__HELP__;
=====================================================================================================
Usage:	Perform KEGG enrichment analysis for noref transcriptome programs.
	perl $0 -diff <..> -ko <..> -g <..> -out <..> -m <..> -padj <..> -c <..> -help

	-diff		different expressed gene list file.
	-ko		kO annotation file produced by annotate.py (in kobas).

	-g		group name, this is the prefix of enrichment result files.
			Defalt is the parameter of -diff.
	-out		output dir. Defalt is the current dir.


	-m		choose statistic method, b is binomial test, c is chi-
                        square test, f is fisher exact test, h is
                        hypergeometric test and x is frequency list, default
                        hypergeometric test.
	-padj		choose false discovery rate (FDR) correction method:
                        QVALUE, BH, BY or none, default QVALUE.
	-c		the cutoff of the gene numbers in a term, default 5.
	-help		show this help message and exit.

Author:	zhangmin
Email:	zhangmin\@novogene.cn
Date:	11/1/2012	ver1.

modified: 10/12/2014 cuijie  #update to kobas2.0-20140801#

=====================================================================================================
__HELP__
}
