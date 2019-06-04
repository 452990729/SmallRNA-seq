#!/bin/perl -w
use strict;
use Getopt::Long;
my ($help, $sample, $type, $adir, $contractId,$code_number);
GetOptions(
	"h|help"	=>\$help,
	"adir=s"	=>\$adir,
	"s|sample=s"	=>\$sample,
	"cont=s"	=>\$contractId,
        "code_num=s"    =>\$code_number,
);
my $usage=<<END;
-----------------------------------------------------------
	perl $0 
	-h|help		help
	-adir=s		Analysis result dir
	-s|sample=s	sample name, splited by ","
	-cont=s		contractId
        -code_num=s     code_number
-----------------------------------------------------------
END
die $usage if ($help or !$adir or !$sample or !$contractId or !$code_number);

my @samples=split/,/,$sample;
my $data_give="$adir/data_give/$code_number\_data_give";
if(!-e "$data_give")
{
	`mkdir -p $data_give/{raw_data,clean_data}`;
}

my $qc="$adir/1.QC/";

## raw
foreach my $i(@samples)
{
	`cp $qc/$i/${i}.raw.fastq.gz  $data_give/raw_data/${i}_raw.fq.gz`;
}	
chdir "$data_give/raw_data";
`md5sum *raw.fq.gz >raw_md5sum.txt`;
chdir "$data_give";
## clean
foreach my $i(@samples)
{
	chdir "$qc/$i/clean_data";
	if (!-e "${i}_clean.fa.gz"){
		`gzip  ${i}_clean_total.fa `;
		`mv ${i}_clean_total.fa.gz ${i}_clean.fa.gz`;
	}
	`mv $qc/$i/clean_data/${i}_clean.fa.gz $data_give/clean_data/${i}_clean.fa.gz`;	
}	
chdir "$data_give/clean_data";
`md5sum *clean.fa.gz >clean_md5sum.txt`; 

chdir "$data_give";
`awk '{print \$1"  raw_data/"\$2}' raw_data/raw_md5sum.txt >> md5sum.txt`;
`awk '{print \$1"  clean_data/"\$2}' clean_data/clean_md5sum.txt >> md5sum.txt`;
`cp /PUBLIC/source/RNA/release_file/验证文件完整性.exe  .`;
`cp /PUBLIC/source/RNA/release_file/数据完整性检验工具使用说明.pdf .`;
#`rm -rf raw_data`;
#`rm -rf clean_data`;

## results
chdir "$data_give";
`mv $adir/$contractId\_sRNA_result/result.py ../../.`;
`mv $adir/$contractId\_sRNA_result/report.py ../../.`;
`rm $adir/$contractId\_sRNA_result/result.log`;
`rm $adir/$contractId\_sRNA_result/report.log`;
`ln -sf $adir/$contractId\_sRNA_result/ `;
`tar -h -czvf $contractId\_sRNA_result.tar.gz  $contractId\_sRNA_result`;
`rm $contractId\_sRNA_result`;
#lib_path.xls
`awk 'BEGIN{FS="\\t";OFS="\/"}{print \$1,\$2}' $adir/mapfile > $adir/data_give/data_path.xls`;
#scripts tar
chdir "$adir";
`nohup sh /PUBLIC/source/RNA/tools/tP &`;
`mv records.tar.gz $adir/data_give`;
#`rm $contractId`;

