#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use FindBin '$Bin';

my ($help, $adir, $outdir, $contID);
GetOptions(
	"h|help" =>\$help,
	"adir=s" =>\$adir,
	"od:s"	=>\$outdir,
	"cont=s" =>\$contID
);
my $usage=<<END;
---------------------------------------------------------------
perl $0 
-h|help	help
-adir=s	analysis dir
-od:s	outputdir, default "/PUBLIC/source/RNA/smallRNA/version3/DATAINFO/yourcontractID"
-cont=s	contractID
---------------------------------------------------------------
END
die $usage if ($help or !$adir or !$contID);
$outdir||="$Bin/../../BioDB/ProjRecord/$contID";
mkdir $outdir if !-e $outdir;
$outdir=abs_path($outdir);
my $outfile="$outdir/$contID\_sRNA_retain_list";
`>$outfile`;

## table
`echo -e "## 1. Raw data basic info" >>$outfile`;
`cat $adir/*QC/*QC_results/results/2RawData_Stat/RawData_Stat.xls >>$outfile`;

`echo -e "\n## 2. Raw data filter info" >>$outfile`;
`cat $adir/*QC/*QC_results/results/3ReadsClassification/clean_process_overview.xls >>$outfile `;

`echo -e "\n## 3. sRNA length filter info" >>$outfile`;
`cat $adir/*QC/*QC_results/results/4length_filter/total_uniq.xls >>$outfile `;

`echo -e "\n## 4. Mapping rate " >>$outfile`;
`cat $adir/*map/*/reference.mapping.stat >>$outfile `;

`echo -e "\n## 5. known miRNA number info " >>$outfile`;
`cat $adir/*known/*/*.known/known_miRNA.map.stat >>$outfile `;

`echo -e "\n## 6. novel miRNA number info " >>$outfile`;
`cat $adir/*novel/*/*.novel/novel_miRNA.map.stat >>$outfile `;

`echo -e "\n## 7.1 sRNA category rc.info " >>$outfile`;
`cat $adir/*Category/category_rc_full.txt >>$outfile `;

`echo -e "\n## 7.2 sRNA category uc.info " >>$outfile`;
`cat $adir/*Category/category_uc_full.txt >>$outfile `;

`echo -e "\n## 8. sRNA diff info " >>$outfile`;
`wc -l $adir/*diff/*/*vs*/*DElist_up.txt | while read num file; do echo -e \$(basename \$file)"\t"\$num; done >$outdir/up`;
`wc -l $adir/*diff/*/*vs*/*DElist_down.txt | while read num file; do echo -e  \$(basename \$file)"\t"\$num; done >$outdir/down`;
`paste $outdir/up $outdir/down >>$outfile`;
#`rm -rf $outdir/up $outdir/down`; 

## image
`cp $adir/*QC/*QC_results/results/4length_filter/*png $outdir`;
`cp $adir/*diff/*/corr_plot/cor_pearson.png $outdir`;

