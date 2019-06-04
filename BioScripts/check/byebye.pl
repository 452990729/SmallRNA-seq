#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
my ($help, $sample, $type, $adir, $organism);
GetOptions(
	"h|help"	=>\$help,
	"adir=s"	=>\$adir,
	"s|sample=s"	=>\$sample,
	"org=s"	=>\$organism,
);
my $usage=<<END;
-----------------------------------------------------------
	perl $0 
	-h|help		help
	-adir=s		Analysis result dir
	-s|sample=s	sample name, splited by ","
	-org=s		organism, only be "plant" or "animal"
-----------------------------------------------------------
END
die $usage if ($help or !$adir or !$sample or !$organism);
my @samples=(split/,/,$sample);

## 1.QC
`rm -r $adir/1.QC/*_vs_*`;  #删除样本vs样本venn 图
for my $sam (@samples){
	!system"rm -r $adir/1.QC/$sam/Category" or print STDERR "$adir/1.QC/$sam/Category is already deleted\n";
	chdir "$adir/1.QC/$sam/clean_data/";
	# 只留下长度筛选后的数据；长度分布图片；clean data;
	if ( ! -e "$sam\_remain_total.fa.gz" ){
		!system"gzip $sam\_remain_total.fa" or die " NOTE: $adir $sam\_remain_total.fa is not exist";
	}
	if ( ! -e "$sam\_remain_uniq.fa.gz" ){
		!system"gzip $sam\_remain_uniq.fa" or die "NOTE: $adir $sam\_remain_uniq.fa is not exist";
	}
	if ( ! -e "$sam\_clean.fa.gz" ){
		system"gzip -c  $sam\_clean_total.fa >$sam\_clean.fa.gz";
	}
	#保留的文件；
	my $qc={
		"$sam\_clean.fa.gz" => 1,
		"$sam\_seq_len_distribution.png" => 1,
		"$sam\_seq_len_distribution.pdf" => 1,
		"$sam\_remain_total.fa.gz" => 1,
		"$sam\_remain_uniq.fa.gz" => 1,
		};
	my $file=[glob('*')];
	&delete($file,$qc);	
}
my $abbr=basename((glob "$adir/2.map/???")[0]);
{## 2.map

#`echo rm -r $adir/2.map/*/tmp >>$adir/byebye_list.sh`;
#	!system"rm -r $adir/2.map/$abbr/tmp" or print "$adir/2.map/$abbr/tmp is already deleted\n";
	chdir "$adir/2.map/$abbr";
	my $map={
		"reference.mapping.stat" => 1,
		"*.map.sh" => 1,
		};
	my $file=[glob'*'];
	&delete($file,$map);

}

{##known and novel
	`rm  $adir/3.known/ref_hairpin.fa`;
	`rm  $adir/3.known/ref_mature.fa`;
	`rm  $adir/3.known/$abbr/ref/*`;
	`rm  $adir/*novel/$abbr/ref/*`;
	`rm -r $adir/*novel/$abbr/$abbr.predict`;
	`rm $adir/3.known/$abbr/$abbr.known/miRBase.mrd`;
	`rm $adir/*novel/$abbr/$abbr.novel/miRBase.mrd`;
	`rm $adir/*novel/$abbr/$abbr.predict/predict.output_permuted.mrd`;
	`rm -r $adir/*novel/$abbr/expression_analyses`;
	`rm -r $adir/*known/$abbr/expression_analyses`;
	`rm -r $adir/*novel/$abbr/pdfs_hsa.novel`;
	`rm -r $adir/*known/$abbr/pdfs_hsa.known`;
}
{
## 3.repeat
my $repeat_dir=(glob"$adir/*repeat")[0];
chdir $repeat_dir;
my @png=glob "*.png";
if(@png>0)  # ref
	{
	#`echo rm -r $adir/5.repeat/*/ >>$adir/byebye_list.sh`;
	my $repeat={
		"repeat.unmap.fas" => 1,
		"repeat.map.fas" => 1,
		"generate_repeat.sh" => 1,
		"$abbr\_rep_map.sh" => 1,
		"$abbr\_rep_map.sh.e" => 1,
		};
	chdir $abbr;
	my $file=[glob'*'];
	&delete($file,$repeat);
	}
else {# noref
	`mv $adir/*.repeat/*/predict/*.gff $adir/5.repeat`;
	`mv $adir/*.repeat/*/predict/DR_TR.fna $adir/5.repeat`;
	`rm -r $adir/*.repeat/$abbr/predict`;
	`rm -r $adir/*.repeat/repeat.unmap.fas`;
	}
}

{#ncRNA
	`rm -r $adir/*.ncRNA/$abbr/input/`;
	chdir "$adir/4.ncRNA/$abbr/output/";
	my $ncRNA={
		"ncRNA.unmap.fas" => 1,
	};
	my $file=[glob'*'];
	&delete($file,$ncRNA);

}

{## 4.target
`rm -r $adir/*target/*/*.fasta` if ($organism=~ /refanimal/);
`cp $adir/*target/*/miRanda/miranda_targets.pairs $adir/*target/*/Common/ ` if ($organism=~ /refanimal/);
`gzip $adir/*target/*/Common/miranda_targets.pairs` if ($organism=~ /refanimal/);
`cp $adir/*target/*/miRanda/miranda_targets_out.fmt $adir/*target/*/Common/ ` if ($organism=~ /refanimal/);
`gzip $adir/*target/*/Common/miranda_targets_out.fmt` if ($organism=~ /refanimal/);
`cp $adir/*target/*/PITA/PITA_pita_results_targets.tab $adir/*target/*/Common/ ` if ($organism=~ /refanimal/);
`gzip $adir/*target/*/Common/PITA_pita_results_targets.tab ` if ($organism=~ /refanimal/);
`cp $adir/*target/*/RNAhybrid/RNAhybrid_miRNA_target_pairs $adir/*target/*/Common/ ` if ($organism=~ /refanimal/);
`gzip $adir/*target/*/Common/RNAhybrid_miRNA_target_pairs ` if ($organism=~ /refanimal/);
`rm -r $adir/*target/*/miRanda/*` if ($organism=~ /refanimal/);
`rm -r $adir/*target/*/PITA/*` if ($organism=~ /refanimal/);
`rm -r $adir/*target/*/RNAhybrid/*` if ($organism=~ /refanimal/);
`gzip $adir/*target/$abbr/*miranda_targets_out.fmt` if ($organism=~ /refanimal/);
`rm -r $adir/*target/*/*miranda_targets_out` if ($organism =~ /refanimal/);
`rm -r $adir/*target/*/*miranda_targets` if ($organism =~ /refanimal/);
`rm -r $adir/*target/*/*_mature.miranda_targets.pairs.annotate` if ($organism =~ /refanimal/);
`rm -r $adir/*target/*/*_mature.miranda_targets.pairs.mg` if ($organism =~ /refanimal/);
`rm -r $adir/*target/*/*_mature.miranda_targets_transcript.pairs` if ($organism =~ /refanimal/);
}


{
## 6.gene
	`rm -r $adir/*gene/$abbr/input`;
	my $gene_dir=(glob"$adir/*gene/$abbr/output")[0];
	print "$gene_dir\n";
	chdir "$gene_dir";
	my $gene={
		"gene.unmap.fas" => "1",
		};
	my $file=[glob'*'];
	print Dumper $file;
	&delete($file,$gene);
}
##edit_family
`rm -r $adir/*edit_family/$abbr/expression_analyses` ;


##novel 过程文件
`rm -r  /PUBLIC/software/RNA/srna-tools-cli/jobs/* ` ;


##enrich
 `rm -r $adir/*enrich/$abbr`;

{## NAT AND TAS
 my $NAT=(glob "$adir/*NAT")[0];
 if ( -e "$NAT" ){
	
	`rm -r $adir/6.NAT/*/*/*.ebwt`;
	my $nat={
		"NAT.unmap.fas" => 1,
		};
	chdir "$adir/6.NAT/$abbr/output/";
	my $file=[glob'*'];
	&delete($file,$nat);
 }
 my $TAS=(glob"$adir/*TAS")[0];
 if (-e $TAS){
	`rm -r $TAS/$abbr/known_TAS/ref`;
	`rm -r $TAS/*/*/*.ebwt`;
	`rm -r $TAS/$abbr/phase.predict/clean_total.fa`;
	chdir "$TAS/$abbr/output";
	my $tas={
		"TAS.unmap.fas" => 1,
		};
	my $file=[glob'*'];
	&delete($file,$tas);
 }


`rm -r $adir/Blast`;
}
sub delete {
	my $a=shift;
	my $b=shift;
	for my $file ( @$a ){
		if ( exists $b->{$file} ){
			next;
		}else{
			!system "rm -r $file" or print STDERR "$adir  $file already deleted\n";
		}
	}
}

