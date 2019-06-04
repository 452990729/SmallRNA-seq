#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Config::Tiny;

sub usage
{
        print STDERR <<USAGE;
==================================================================================
Description  run DEGseq and plotting
Version:1.0
2013-01-04  jiangxiaoxue\@novogene.cn

Options
<Required>
	-r <s>  :	1.readcount,2.readcount
			the readcount matrix from known miRNA analysis and novel miRNA analysis or others,seperated by ","
	-s <s>  :       samplenames , sep=","(such as sample1,sample2,sample3...)
<Optional>[default, only suitable for two samples]
	-group <s> :	group that the samples belong, default=(-s)
	e.g. sample1,sample2:sample3(sample2 and sample3 are the same group)
	-groupname <s>:	group names, group1,group2,... [split by ",", default=(-s)]
	-g <s>  :	how to compare the groups, e.g. 1:2,1:3; default 1:2(test:control)
			[if more than two sample, please makesure defined this parameters]
	-venn <s> :	how to plot venn ; eg : inorder to get the venn of comparision 1:2, 1:3and 2:3, and get another venn of diffgeneIDs of comparision 1:4 and 2:3, you need to input 1:2_1:3_2:3,1:4_2:3
	-o <s>	:	the output project, default ./DE.out
	-p1 <f>	:	DESeq: the threshold of adjusted pvalue for diff sRNA, default=0.05
	-p2 <f>	:	DEGseq: the threshold of adjusted pvalue for diff sRNA, default=0.01
	-f1 <f>	:	DESeq: the threshold of |Foldchange| for diff sRNA, default=2
	-f2 <f>	:	DEGseq: the threshold of |Foldchange| for diff sRNA, default=2
	-h|?|help:	Show this help
===============================================================================
USAGE
}
my($rc, $samples, $group, $groupname, $compare, $venn, $dir, $p1val, $Fc1, $p2val, $Fc2, $help);
GetOptions(
        "h|?|help"=>\$help,
        "r=s"=>\$rc,
	"s=s"=>\$samples,
	"group=s"=>\$group,
	"groupname=s"=>\$groupname,
	"g=s"=>\$compare,
	"venn=s"=>\$venn,
	"o=s"=>\$dir,
	"p1=f"=>\$p1val,
	"p2=f"=>\$p2val,
	"f1=f"=>\$Fc1,
	"f2=f"=>\$Fc2,
);

if((!defined($rc))||(!defined($samples))|| (defined($help))){
        &usage;
        exit 0;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};

$p1val ||=0.05;
$p2val ||=0.01;
$Fc1 ||=2;
$Fc2 ||=2;
my @sample=split /,/,$samples;

my @group_compose;
my @group_name;
if(!defined($group)||!defined($groupname)){
	if(!defined($group)){$group=$samples;}
	if(!defined($groupname)){$groupname=$samples;}
	if(scalar(@sample)!=2 && !defined($compare)){&usage;exit 0;}
}
@group_compose=split /,/,$group;
@group_name=split /,/,$groupname;
if(scalar(@group_compose)!=scalar(@group_name)){&usage;die "please makesure the corresponding number between \"group composition\" and groupname\n";}

$compare ||= "1:2";
my @com=split /,/,$compare;

my $pwd = `pwd`;
$pwd =~ s/\s+//;
$dir ||= "$pwd/DE.out";
if(!(-e $dir)){
	`mkdir $dir`;
}
$dir=&get_ful_path($dir);

#merged these readcount file, and get tpm
`$perlExec $Bin/sRNA_merge_readcount.pl $rc $dir/merged`;

open(OUT1,">$dir/group.txt");
print OUT1 "sample\tgroup\n";
my %group_compose_detail;
foreach my $i(0..$#group_compose){
	my @tmp=split(":",$group_compose[$i]);
	foreach my $j(@tmp){
		print OUT1 "$j\t$group_name[$i]\n";
		$group_compose_detail{$group_name[$i]}{$j}=1;
	}
}
close(OUT1);

open(OUT2,">$dir/compare.txt");
foreach my $i(@com){
	my @tmp=split(":",$i);
	my $outdir=$group_name[$tmp[0]-1]."vs".$group_name[$tmp[1]-1];
	print OUT2 $outdir,"\t",$group_name[$tmp[0]-1],"\t",$group_name[$tmp[1]-1],"\n";
	$outdir=$dir."/".$outdir;
	if($group_compose[$tmp[0]-1]=~/:/||$group_compose[$tmp[1]-1]=~/:/){
		print "$perlExec $Bin/DESeq2_sRNA.v3.pl -i $dir/meanscount.txt -a $group_compose[$tmp[0]-1] -b $group_compose[$tmp[1]-1] -n1 $group_name[$tmp[0]-1] -n2 $group_name[$tmp[1]-1] -op $outdir -p $p1val\n";
	}else{
		print "$perlExec $Bin/DEGseq.R.v3.pl -i $dir/meanscount.txt -r $dir/meanstpm.txt -a $group_compose[$tmp[0]-1] -b $group_compose[$tmp[1]-1] -o $outdir -p $p2val -f $Fc2\n";
	}
}
print "\nsh diffsum.sh\n";
close(OUT2);

print "\n### Density, boxplot, corr plot, cluster (Hcluster,Kmeans,SOM)\n";
print "$perlExec $Bin/plot.R.v2.pl -rc $dir/merged.readcount -mean_rc $dir/meanscount.txt -tpm $dir/merged.tpm -mean_tpm $dir/meanstpm.txt -compare $dir/compare.txt -output $dir -cluster $samples\n\n";
if (@group_name == 2 and not defined $venn){
	print "$perlExec $Bin/samples_coExpression_venn2_replicates.pl -venn $group_name[0]:$group_name[1] -fpkm $dir/meanstpm.txt -out $dir/venn\n";
}#edit in 2015/3/4 
if(defined($venn)){
	print "$perlExec $Bin/plot_only_venn.pl -in $dir -r $dir/merged.tpm -compare $dir/compare.txt -groupnames $groupname -venn $venn -output $dir/venn\n";
}

=head
print "\n### miRNA correlation heatmap\n";
print "cd $dir/corr_plot\n";
print "perl $Bin/bin/miRNA_corheatmap.pl -tpm $dir/merged.tpm\n";
print "mv miRNA_TPM_corheat.R ../";
=cut


open(IN1,"$dir/merged.readcount");
my $header1=<IN1>;
chomp($header1);
my @header1=(split/\t/,$header1); shift @header1; # 2013-10-22 for check yangzie
my %count;
my %count_sample;  # 2013-10-22 for check yangzie
my ($types1,@sample_title1)=split("\t",$header1);
while(<IN1>){
	chomp;
	my @tmp=split /\t/;
	for my $i(1..$#tmp){
		$count{$tmp[0]}{$sample_title1[$i-1]}=$tmp[$i];
	}
	$count_sample{$tmp[0]}=join("\t",@tmp[1..$#tmp]); # 2013-10-22 for check yangzie
}
close(IN1);
open(IN2,"$dir/merged.tpm");
my $header2=<IN2>;
chomp($header2);
my @header2=(split/\t/,$header2); shift @header2; # 2013-10-22 for check yangzie
my %tpm;  
my %tpm_sample; #  2013-10-22 for check yangzie
my ($types2,@sample_title2)=split("\t",$header2);
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	for my $i(1..$#tmp){
		$tpm{$tmp[0]}{$sample_title2[$i-1]}=$tmp[$i];
	}
	$tpm_sample{$tmp[0]}=join("\t",@tmp[1..$#tmp]); # 2013-10-22 for check yangzie
}
close(IN2);

foreach my $i(keys %group_compose_detail){
	my $number=scalar(keys $group_compose_detail{$i});
	if($number!=1){
		foreach my $j(keys $group_compose_detail{$i}){
			foreach my $m(keys %count){
				$count{$m}{$i}+=$count{$m}{$j}/$number;
				$tpm{$m}{$i}+=$tpm{$m}{$j}/$number;
			}
		}
	}
}
open(OUT1,">$dir/meanscount.txt");
open(OUT2,">$dir/meanstpm.txt");
=head
print OUT1 "sRNA\t",join("\t",@group_name),"\t",join("\t",@header1),"\n"; # 2013-10-22 for check yangzie
print OUT2 "sRNA\t",join("\t",@group_name),"\t",join("\t",@header2),"\n"; # 2013-10-22 for check yangzie
foreach my $i(sort {$a cmp $b} keys %count){
	print OUT1 "$i";
	print OUT2 "$i";
	foreach my $j(@group_name){
		print OUT1 "\t",$count{$i}{$j};
		print OUT2 "\t",$tpm{$i}{$j};
	}
	print OUT1 "\t$count_sample{$i}\n"; # 2013-10-22 for check yangzie
	print OUT2 "\t$tpm_sample{$i}\n"; # 2013-10-22 for check yangzie
}
close(OUT1);
close(OUT2);
=cut 
my @name;   # No need %count_sample!  2013-10-24 modify
for (sort keys %count)
{
    @name=sort keys $count{$_};
}
print OUT1 "sRNA\t",join("\t",@name),"\n";
print OUT2 "sRNA\t",join("\t",@name),"\n";
for my $sRNA(sort keys %count)
{
    print OUT1 "$sRNA";
    print OUT2 "$sRNA";
    for my $mem(@name)
    {   
        print OUT1 "\t$count{$sRNA}{$mem}";
        print OUT2 "\t$tpm{$sRNA}{$mem}";
    }   
    print OUT1 "\n";
    print OUT2 "\n";
}

close(OUT1);
close(OUT2);

sub get_ful_path{
	my $in=shift;
	if($in !~ /^\//){
		my $t=`pwd`;
		chomp($t);
		return "$t/$in";
	}else{
		return $in;
        }
}

