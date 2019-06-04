use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use File::Basename;

sub usage{
print STDERR <<USAGE;
=========================================================================
Usage: perl $0 [options]
example:perl $0 -n sp1,sp2 -dir QC_dir -min 18 -max 40  -cont NH*
options:
[mandatory parameters]
	-n sp1,sp2      List of sample name(3 character) corresponding to reads,seperated by ","
	-cont			contratID
	-dir	<s>		QC directory of this project
	-min	<i>		Minimum length of sRNA [18]
	-max	<i>		Maximum length of sRNA [40]
	-mapfile <s>	mapfile 3column: "Path,Library\tSampleName"
[optional parameters]
	-o	out_directory[default="projpath/QC_results"]
	-h|?|help   :  Show this help
=========================================================================
USAGE
}
my ($samples,$dir,$out,$min,$max,$help,$contractId, $mapfile);
GetOptions(
	"h|?|help"=>\$help,
	"cont=s"=>\$contractId,
	"n=s"=>\$samples,
	"dir=s"=>\$dir,
	"o=s"=>\$out,
	"min=i"=>\$min,
	"max=i"=>\$max,
	"mapfile=s" =>\$mapfile,
);

$min = 18  if(!defined($min));
$max = 40  if(!defined($max));
if(!defined($samples) || !defined($dir) || !defined($mapfile) || defined($help)){
	&usage;
	exit 0;
}
my @sample=(split /,/, $samples);

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
my $projpath = dirname($dir);
$out||="$projpath/$contractId\_QC_results";
if(!(-e $out)){ `mkdir -p $out`; }
`cp -r $Bin/src $out`;

if(!(-e "$out/results")){ `mkdir -p $out/results`; }
my ($head,$string);

my %hash;
open M, $mapfile or die $!;
while(<M>)
{
	chomp;
	my($path,$library, $smpName)=split/\t/;
	$hash{$smpName}=$library;
}
close M;

## 1. Raw data Errorate
if(!(-e "$out/results/1RawData_ErrorRate")){
	`mkdir -p $out/results/1RawData_ErrorRate`;
}
my $errorpng;
for my $i(@sample)
{
    `cp $dir/$i/clean_data/$i.clean_Error.png $out/results/1RawData_ErrorRate/$i\_error_rate_distribution.png`;
    `convert -resize 600 $dir/$i/clean_data/$i.clean_Error.png $out/src/images/$i.error_rate_distribution.png`;
    `convert -resize 90 $dir/$i/clean_data/$i.clean_Error.png $out/src/images/$i.error_rate_distribution.JPEG`;
    $errorpng.=pngerror($i);
}

## 2. Raw data basic info 
if(!(-e "$out/results/2RawData_Stat")){
	`mkdir -p $out/results/2RawData_Stat`;
}
$string="";
$head="Library\tSample\tReads\tBases\tError rate\tQ20\tQ30\tGC content\n";
my $tab_2_1=headtb2htmltb($head);
#$html.=headtb2htmltb($head);
foreach my $i(@sample){
	my $tmp=`tail -1 $dir/$i/clean_data/$i\_raw.stat`;
	my $smpName=$hash{(split/\t/,$tmp)[0]};
	$tmp=$smpName."\t".$tmp;
	#$html.=tb2htmltb($tmp);
	$tab_2_1.=tb2htmltb($tmp);
	$string.="$tmp\n";
}
`perl -e 'print "$head$string"' >$out/results/2RawData_Stat/RawData_Stat.xls`;

## 3. Raw data QC filter
if(!(-e "$out/results/3ReadsClassification")){
	`mkdir -p $out/results/3ReadsClassification`;
}
$string="";
$head="Library\tSample\ttotal_reads\tN% > 10%\tlow quality\t5_adapter_contamine\t3_adapter_null or insert_null\twith ployA/T/G/C\tclean reads\n";
my $tab_3_1.=headtb2htmltb($head);
foreach my $i(@sample){
	my $tmp=$i;
	$tmp=$hash{$i}."\t".$i;
	open(IN,"$dir/$i/clean_data/$i\_clean_process_overview.txt");
	my $recording=0;
	while(<IN>){
		chomp;
		$recording++;
		if($recording!=6 && ($recording<9)){
			my @temp=split /\t/;
			$temp[2] =~ s/[()% ]//g;
			$temp[2] = sprintf("%.2f",$temp[2]);
			$tmp.="\t$temp[1] ($temp[2]%)";
		}
	}
	$tmp.="\n";
	$tab_3_1.=tb2htmltb($tmp);
	$string.=$tmp;
}
`perl -e 'print "$head$string"' >$out/results/3ReadsClassification/clean_process_overview.xls`;

## 4. Length distribution stat and plot 
if(!(-e "$out/results/4length_filter")){
	`mkdir -p $out/results/4length_filter`;
}
$string="";
$head="Library\tSample\tTotal reads\tTotal bases (bp)\tUniq reads\tUniq bases (bp)\n";
my $tab_4_1.=headtb2htmltb($head);
foreach my $i(@sample){
	my $tmp=$i;
	$tmp=$hash{$i}."\t".$i;
	open(IN,"$dir/$i/clean_data/$i\_remain_total_uniq.txt");
	while(<IN>){
		chomp;
		my @temp=split /\t/;
		$tmp.="\t$temp[1]";
	}
	$tmp.="\n";
	#$html.=tb2htmltb($tmp);
	$tab_4_1.=tb2htmltb($tmp);
	$string.=$tmp;
}
`perl -e 'print "$head$string"' >$out/results/4length_filter/total_uniq.xls`;

my $fig_4_2;
foreach my $i(@sample){
	`cp $dir/$i/clean_data/$i\_len.stat $out/results/4length_filter/$i\_seq_len_distribution.txt`;
	`cp $dir/$i/clean_data/$i\_seq_len_distribution.png $out/results/4length_filter/$i\_seq_len_distribution.png`;
	`cp $dir/$i/clean_data/$i\_seq_len_distribution.pdf $out/results/4length_filter/$i\_seq_len_distribution.pdf`;
	$fig_4_2.="<img class=\"w45\" src=\"results/4length_filter/$i\_seq_len_distribution.png\" />";
}

my $lengthpng;
for my $i(@sample)
{
    `convert -resize 600 $dir/$i/clean_data/$i\_seq_len_distribution.png $out/src/images/$i\_seq_len_distribution.png`; 
    `convert -resize 90 $dir/$i/clean_data/$i\_seq_len_distribution.png $out/src/images/$i\_seq_len_distribution.JPEG`; 
    $lengthpng.=pnglength($i);
}

## 5. mapping rate
my $map_dir=$dir;
$map_dir=~ s/1\.QC//g;
$map_dir =glob "${map_dir}2.map/reference.mapping.stat";
open MAP,$map_dir;
$head=<MAP>;
chomp($head); 
my $tab_5_1 =headtb2htmltb($head);
while (<MAP>) {
	chomp;
	$tab_5_1 .=tb2htmltb($_);
}
close MAP;
## 6. Rfam aligning
if(!(-e "$out/results/5category/")){
	`mkdir -p $out/results/5category`;
}
$string="";
my $tab_6_1;
foreach my $i(@sample){
	my $tmp=$i;
	$tmp=$hash{$i}."\t".$i;
	my $id=0;
	open(IN,"$dir/$i/Category/$i\_category.log");
	while(<IN>){
		chomp;
		$id++;
		if($id>5 && $id<12){
			my @temp=split /\t/;
			$tmp.="\t$temp[0]($temp[2],$temp[3])";
		}
	}
	$tmp.="\n";
	#$html.=tb2htmltb($tmp);
	$tab_6_1.=tb2htmltb($tmp);
	$string.=$tmp;
}
`perl -e 'print "$string"' >$out/results/5category/sRNA_Rfam_category_top5.xls`;


# ************************************ html  ***********************************************
my $html=<<END;
<!DOCTYPE html PUBLIC "-//W3C//DTD Xhtml 1.0 Transitional//EN">
<html>
<head>
<title>诺禾致源Small RNA生物信息分析报告</title>
<META NAME="Author" CONTENT="wangshaobin\@novogene.cn"> 
<META NAME="Version" CONTENT="20140114">
<meta http-equiv="content-type" content="text/html;">
<link rel="stylesheet" type="text/css" href="src/css/text.css">
<link rel="StyleSheet" href="src/js/tree/tree.css" type="text/css">
<link rel="stylesheet" type="text/css" href="src/js/fancybox/jquery.fancybox-1.3.4.css" media="screen" />
<link rel="stylesheet" href="src/css/style.css" />
<script src="src/js/jquery-1.4.2.min.js" type="text/javascript"></script>
<script src="src/js/scrollTop.js" type="text/javascript"></script>
<script src="src/js/common.js" type="text/javascript"></script>
<script src="src/js/jquery.albumSlider.min.js" type="text/javascript"></script>
<script type="text/javascript" src="src/js/fancybox/jquery.mousewheel-3.0.4.pack.js"></script>
<script type="text/javascript" src="src/js/fancybox/jquery.fancybox-1.3.4.pack.js"></script>
<script type="text/javascript" src="src/js/tree/tree.js"></script>

<script type="text/javascript">
	\$(document).ready(function() {
		/*
		*   Examples - images
		*/

		\$("a#example1").fancybox();

		\$("a#example2").fancybox({
			'overlayShow'	: false,
			'transitionIn'	: 'elastic',
			'transitionOut'	: 'elastic'
		});
		\$("a#example3").fancybox({
			'transitionIn'	: 'none',
			'transitionOut'	: 'none'	
		});

		\$("a#example4").fancybox({
			'opacity'		: true,
			'overlayShow'	: false,
			'transitionIn'	: 'elastic',
			'transitionOut'	: 'none'
		});

		\$("a#example5").fancybox();

		\$("a#example6").fancybox({
			'titlePosition'		: 'outside',
			'overlayColor'		: '#000',
			'overlayOpacity'	: 0.9
		});

		\$("a#example7").fancybox({
			'titlePosition'	: 'inside'
		});

		\$("a#example8").fancybox({
			'titlePosition'	: 'over'
		});

		\$("a[rel=example_group]").fancybox({
			'transitionIn'		: 'none',
			'transitionOut'		: 'none',
			'titlePosition' 	: 'over',
			'titleFormat'		: function(title, currentArray, currentIndex, currentOpts) {
				return '<span id="fancybox-title-over">Image ' + (currentIndex + 1) + ' / ' + currentArray.length + (title.length ? ' &nbsp; ' + title : '') + '</span>';
			}
		});

		/*
		*   Examples - various
		*/

		\$("#various1").fancybox({
			'titlePosition'		: 'inside',
			'transitionIn'		: 'none',
			'transitionOut'		: 'none'
		});

		\$("#various2").fancybox();

		\$("#various3").fancybox({
			'width'				: '75%',
			'height'			: '75%',
			'autoScale'			: false,
			'transitionIn'		: 'none',
			'transitionOut'		: 'none',
			'type'				: 'iframe'
		});

		\$("#various4").fancybox({
			'padding'			: 0,
			'autoScale'			: false,
			'transitionIn'		: 'none',
			'transitionOut'		: 'none'
		});
	});
</script>

<style media="print">
.noprint {DISPLAY: none;}
</style>

<div class="noprint">
<div style="display: block;" id="goTopBtn">
<a class="backtotop" title="回顶部"><img src="src/images/goTop.jpg" width="30" height="30" class="back-tip"/></a>
</div>
</div>

<script type="text/javascript"> 
function displaySubMenu(li) { 
var subMenu = li.getElementsByTagName("ul")[0]; 
subMenu.style.display = "block"; 
} 
function hideSubMenu(li) { 
var subMenu = li.getElementsByTagName("ul")[0]; 
subMenu.style.display = "none"; 
} 
</script> 

	<script type="text/javascript">
		var Tree = new Array;
		// nodeId | parentNodeId | nodeName | nodeUrl
	<{jtree}>	
	</script>
<script type="text/javascript">
\$(function(){
	//纵向，默认，移动间隔2
    \$('div.albumSlider').albumSlider();
    //横向设置
    \$('div.albumSlider-h').albumSlider({direction:'h',step:3});
});
</script>	
</head>

<body>
<!---------------------------------------- 目录 ---------------------------------------------->
<div id="page">
<p><a name="home"><img class="normal" src="src/images/logo.png" /></a></p>

<h1>Small RNA QC分析报告</h1>
<p class="paragraph">
	<ul class="alt">
		<li><a href="#建库测序实验流程">1	测序错误率分布检查</a></li>
		<li><a href="#原始测序数据质量情况汇总">2	原始测序数据质量情况汇总</a></li>
		<li><a href="#测序数据过滤">3	测序数据过滤</a></li>
		<li><a href="#sRNA长度筛选">4	sRNA长度筛选</a></li>
		<li><a href="#参考序列比对情况汇总">5 参考序列比对情况汇总</a></li>
		<li><a href="#sRNA与Rfam数据库比对">6	sRNA与Rfam数据库比对</a></li>
	</ul>
	<br />
</p>
</div>

<!----------------------------------------- 1.测序错误率分布检查 -------------------------------------->
<div id="page">
<p class="head">
	<a href="#home" title = "返回首页"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="测序错误率分布检查">北京诺禾致源生物信息科技有限公司</a> 
	<hr/> 
	<br/>
</p>
<h3>1	测序错误率分布检查</h3>

<div class="albumSlider">
  <div class="fullview"> <img src="src/images/$sample[0].error_rate_distribution.png" /></div>
  <div class="slider">
    <div class="button movebackward" title="向上滚动"></div>
    <div class="imglistwrap"><ul class="imglist">
    $errorpng
    </ul></div>
    <div class="button moveforward" title="向下滚动"></div>
  </div>
</div>

<!----------------------------------------- 2.原始测序数据质量情况汇总 -------------------------------------->
<div id="page">
<p class="head">
	<a href="#home" title = "返回首页"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="原始测序数据质量情况汇总">北京诺禾致源生物信息科技有限公司</a> 
	<hr/>
	<br/>
</p>
<h3>2	原始测序数据质量情况汇总</h3>
<table>
$tab_2_1
</table></p>
<br />
</div>

<!----------------------------------------- 3.测序数据过滤 -------------------------------------->
<p class="head">
	<a href="#home" title = "返回首页"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="测序数据过滤">北京诺禾致源生物信息科技有限公司</a> 
	<hr/> 
	<br/>
</p>
<h3>3 测序数据过滤</h3>
<table>
$tab_3_1
</table></p>
<p class="tremark">注：<br />
(1) Sample：样品id。<br />
(2) total_reads：统计原始序列数据条数。<br />
(3) N% > 10%：因 N含量超过10%，过滤掉的reads数及其占总raw reads数的比例。<br />
(4) low quality：因低质量，过滤掉的reads数及其占总raw reads数的比例。<br />
(5) 5_adapter_contamine：因含有5’接头，过滤掉的reads数及其占总raw reads数的比例。<br />
(6) 3_adapter_null or insert_null：因没有3’接头或没有插入片段，过滤掉的reads数及其占总raw reads数的比例。<br />
(7) with ployA/T/G/C：因含有ployA/T/G/C，过滤掉的reads数及其占总raw reads数的比例。<br />
(8) clean reads：最终得到的clean reads数及其占总raw reads数的比例。</p>
<br/>
</div>

<!----------------------------------------- 4.sRNA长度筛选 -------------------------------------->
<p class="head">
	<a href="#home" title = "返回首页"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="sRNA长度筛选">北京诺禾致源生物信息科技有限公司</a> 
	<hr/> 
	<br/>
</p>
<h3>4	sRNA长度筛选</h3>
<p class="paragraph">对各样品的clean reads，筛选一定长度的sRNA来进行后续分析。以下是对这些小RNA（sRNA）的种类（用uniq表示）及数量（用total表示）<b>【见表5.1】</b>，以及长度分布统计<b>【见图5.2】</b>。一般来说，小RNA的长度区间为18~40nt，长度分布的峰能帮助我们判断小RNA的种类，如miRNA集中在21~22nt，siRNA集中在24nt等。</p>
<p class="name">表5.1 sRNA种类和数量情况一览表 </p>
<table>
$tab_4_1
</table></p>
<p class="tremark">注：<br />
(1) Sample：样品id。<br />
(2) Total reads：sRNA的总数。<br />
(3) Total bases (bp)：sRNA的总长度。<br />
(4) Uniq reads：clean sRNA的种类。<br />
(5) Uniq bases (bp)：各种sRNA的总长度。</p>
<p class="center">

<div class="albumSlider">
  <div class="fullview"> <img src="src/images/$sample[0]\_seq_len_distribution.png" /></div>
  <div class="slider">
    <div class="button movebackward" title="向上滚动"></div>
    <div class="imglistwrap"><ul class="imglist">
    $lengthpng
    </ul></div>
    <div class="button moveforward" title="向下滚动"></div>
  </div>
</div>

</p>
</div>
<!----------------------------------------- 5.参考序列比对分析 -------------------------------------->
<div id="page">
<p class="head">
	<a href="#home" title = "返回首页"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="参考序列比对情况汇总">北京诺禾致源生物信息科技有限公司</a> 
	<hr/>
	<br/>
</p>
<h3>5	参考序列比对情况汇总</h3>
<table>
$tab_5_1
</table></p>
<br />
<p class="tremark">注：<br />
      (1) Sample：样品id。<br />
      (2) Total sRNA：经分析2.4“sRNA长度筛选”后，各个样本所得到的总reads数。<br />
      (3) Mapped sRNA：该样本reads中能mapped到参考序列的reads数及所占百分比。<br />
      (4) “+” Mapped sRNA：该样本reads中能mapped到参考序列方向相同链的reads数及所占百分比。<br />
      (5) “-” Mapped sRNA：该样本reads中能mapped到参考序列方向相反链的reads数及所占百分比。 </p>
</div>
<!----------------------------------------- 6.sRNA与Rfam数据库比对 -------------------------------------->
<div id="page">
<p class="head">
	<a href="#home" title = "返回首页"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="sRNA与Rfam数据库比对">北京诺禾致源生物信息科技有限公司</a> 
	<hr/>
	<br/>
</p>
<h3>6	sRNA与Rfam数据库比对</h3>

<p class="paragraph">将长度筛选后的sRNA与Rfam数据库（nocoding）进行比对，看样本生物学重复之间的sRNA在Rfam数据库上的物种分布情况。若物种分布差异较大，则结合上面的长度分布图，如果差异也较大，则样品QC不达标。</p>
<br />
<table>
$tab_6_1
</table>
</p>
<p class="tremark">注：格式（species,spe_mapped_reads/total_reads,spe_mapped_reads/mapped_reads）（物种，比对上该物种的reads数占总reads的百分比，比对上该物种的reads数占mapped到Rfam整个数据库的reads的百分比）<br />
</p>
<br/>
</div>

</body>
</html>
END

open WEB,">$out/sQCReport.html";
#print WEB encode("utf-8",decode("utf-8",$html));
print WEB $html;
close(WEB);

`chmod a+x $out/sQCReport.html`;

sub headtb2htmltb{
    my $line=shift;
    chomp($line);
    $line = "<tr\> <th>".$line."</th> </tr>\n";
    $line =~ s#\t#</th> <th>#g;
    return $line;
}

sub tb2htmltb{
	my $line=shift;
	chomp($line);
	$line =~ s#\n#</td> </tr>\n<tr> <td>#g;
	$line = "<tr> <td>".$line."</td> </tr>\n";
	$line =~ s#\t#</td> <td>#g;
	return $line;
}

sub pngerror{
    my $line=shift;
    chomp($line);
    $line = "\n\t<li> <a id=\"example2\" href=\'src\/images\/".$line."\.error\_rate\_distribution\.png\'> <img src=\'src\/images\/".$line."\.error\_rate\_distribution\.JPEG\'\/> <\/a> <\/li>";
    return $line;
}

sub pnglength{
    my $line=shift;
    chomp($line);
    $line = "\n\t<li> <a id=\"example2\" href=\'src\/images\/".$line."\_seq_len_distribution.png\'> <img src=\'src\/images\/".$line."\_seq_len_distribution.JPEG\'\/> <\/a> <\/li>";
    return $line;
}

