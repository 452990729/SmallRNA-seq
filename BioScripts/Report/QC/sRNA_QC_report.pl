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
<title>ŵ����ԴSmall RNA������Ϣ��������</title>
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
<a class="backtotop" title="�ض���"><img src="src/images/goTop.jpg" width="30" height="30" class="back-tip"/></a>
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
	//����Ĭ�ϣ��ƶ����2
    \$('div.albumSlider').albumSlider();
    //��������
    \$('div.albumSlider-h').albumSlider({direction:'h',step:3});
});
</script>	
</head>

<body>
<!---------------------------------------- Ŀ¼ ---------------------------------------------->
<div id="page">
<p><a name="home"><img class="normal" src="src/images/logo.png" /></a></p>

<h1>Small RNA QC��������</h1>
<p class="paragraph">
	<ul class="alt">
		<li><a href="#�������ʵ������">1	��������ʷֲ����</a></li>
		<li><a href="#ԭʼ�������������������">2	ԭʼ�������������������</a></li>
		<li><a href="#�������ݹ���">3	�������ݹ���</a></li>
		<li><a href="#sRNA����ɸѡ">4	sRNA����ɸѡ</a></li>
		<li><a href="#�ο����бȶ��������">5 �ο����бȶ��������</a></li>
		<li><a href="#sRNA��Rfam���ݿ�ȶ�">6	sRNA��Rfam���ݿ�ȶ�</a></li>
	</ul>
	<br />
</p>
</div>

<!----------------------------------------- 1.��������ʷֲ���� -------------------------------------->
<div id="page">
<p class="head">
	<a href="#home" title = "������ҳ"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="��������ʷֲ����">����ŵ����Դ������Ϣ�Ƽ����޹�˾</a> 
	<hr/> 
	<br/>
</p>
<h3>1	��������ʷֲ����</h3>

<div class="albumSlider">
  <div class="fullview"> <img src="src/images/$sample[0].error_rate_distribution.png" /></div>
  <div class="slider">
    <div class="button movebackward" title="���Ϲ���"></div>
    <div class="imglistwrap"><ul class="imglist">
    $errorpng
    </ul></div>
    <div class="button moveforward" title="���¹���"></div>
  </div>
</div>

<!----------------------------------------- 2.ԭʼ������������������� -------------------------------------->
<div id="page">
<p class="head">
	<a href="#home" title = "������ҳ"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="ԭʼ�������������������">����ŵ����Դ������Ϣ�Ƽ����޹�˾</a> 
	<hr/>
	<br/>
</p>
<h3>2	ԭʼ�������������������</h3>
<table>
$tab_2_1
</table></p>
<br />
</div>

<!----------------------------------------- 3.�������ݹ��� -------------------------------------->
<p class="head">
	<a href="#home" title = "������ҳ"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="�������ݹ���">����ŵ����Դ������Ϣ�Ƽ����޹�˾</a> 
	<hr/> 
	<br/>
</p>
<h3>3 �������ݹ���</h3>
<table>
$tab_3_1
</table></p>
<p class="tremark">ע��<br />
(1) Sample����Ʒid��<br />
(2) total_reads��ͳ��ԭʼ��������������<br />
(3) N% > 10%���� N��������10%�����˵���reads������ռ��raw reads���ı�����<br />
(4) low quality��������������˵���reads������ռ��raw reads���ı�����<br />
(5) 5_adapter_contamine������5����ͷ�����˵���reads������ռ��raw reads���ı�����<br />
(6) 3_adapter_null or insert_null����û��3����ͷ��û�в���Ƭ�Σ����˵���reads������ռ��raw reads���ı�����<br />
(7) with ployA/T/G/C������ployA/T/G/C�����˵���reads������ռ��raw reads���ı�����<br />
(8) clean reads�����յõ���clean reads������ռ��raw reads���ı�����</p>
<br/>
</div>

<!----------------------------------------- 4.sRNA����ɸѡ -------------------------------------->
<p class="head">
	<a href="#home" title = "������ҳ"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="sRNA����ɸѡ">����ŵ����Դ������Ϣ�Ƽ����޹�˾</a> 
	<hr/> 
	<br/>
</p>
<h3>4	sRNA����ɸѡ</h3>
<p class="paragraph">�Ը���Ʒ��clean reads��ɸѡһ�����ȵ�sRNA�����к��������������Ƕ���ЩСRNA��sRNA�������ࣨ��uniq��ʾ������������total��ʾ��<b>������5.1��</b>���Լ����ȷֲ�ͳ��<b>����ͼ5.2��</b>��һ����˵��СRNA�ĳ�������Ϊ18~40nt�����ȷֲ��ķ��ܰ��������ж�СRNA�����࣬��miRNA������21~22nt��siRNA������24nt�ȡ�</p>
<p class="name">��5.1 sRNA������������һ���� </p>
<table>
$tab_4_1
</table></p>
<p class="tremark">ע��<br />
(1) Sample����Ʒid��<br />
(2) Total reads��sRNA��������<br />
(3) Total bases (bp)��sRNA���ܳ��ȡ�<br />
(4) Uniq reads��clean sRNA�����ࡣ<br />
(5) Uniq bases (bp)������sRNA���ܳ��ȡ�</p>
<p class="center">

<div class="albumSlider">
  <div class="fullview"> <img src="src/images/$sample[0]\_seq_len_distribution.png" /></div>
  <div class="slider">
    <div class="button movebackward" title="���Ϲ���"></div>
    <div class="imglistwrap"><ul class="imglist">
    $lengthpng
    </ul></div>
    <div class="button moveforward" title="���¹���"></div>
  </div>
</div>

</p>
</div>
<!----------------------------------------- 5.�ο����бȶԷ��� -------------------------------------->
<div id="page">
<p class="head">
	<a href="#home" title = "������ҳ"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="�ο����бȶ��������">����ŵ����Դ������Ϣ�Ƽ����޹�˾</a> 
	<hr/>
	<br/>
</p>
<h3>5	�ο����бȶ��������</h3>
<table>
$tab_5_1
</table></p>
<br />
<p class="tremark">ע��<br />
      (1) Sample����Ʒid��<br />
      (2) Total sRNA��������2.4��sRNA����ɸѡ���󣬸����������õ�����reads����<br />
      (3) Mapped sRNA��������reads����mapped���ο����е�reads������ռ�ٷֱȡ�<br />
      (4) ��+�� Mapped sRNA��������reads����mapped���ο����з�����ͬ����reads������ռ�ٷֱȡ�<br />
      (5) ��-�� Mapped sRNA��������reads����mapped���ο����з����෴����reads������ռ�ٷֱȡ� </p>
</div>
<!----------------------------------------- 6.sRNA��Rfam���ݿ�ȶ� -------------------------------------->
<div id="page">
<p class="head">
	<a href="#home" title = "������ҳ"><img class="logo" align="left" src="src/images/logo.png" /></a>
	<a name="sRNA��Rfam���ݿ�ȶ�">����ŵ����Դ������Ϣ�Ƽ����޹�˾</a> 
	<hr/>
	<br/>
</p>
<h3>6	sRNA��Rfam���ݿ�ȶ�</h3>

<p class="paragraph">������ɸѡ���sRNA��Rfam���ݿ⣨nocoding�����бȶԣ�����������ѧ�ظ�֮���sRNA��Rfam���ݿ��ϵ����ֲַ�����������ֲַ�����ϴ���������ĳ��ȷֲ�ͼ���������Ҳ�ϴ�����ƷQC����ꡣ</p>
<br />
<table>
$tab_6_1
</table>
</p>
<p class="tremark">ע����ʽ��species,spe_mapped_reads/total_reads,spe_mapped_reads/mapped_reads�������֣��ȶ��ϸ����ֵ�reads��ռ��reads�İٷֱȣ��ȶ��ϸ����ֵ�reads��ռmapped��Rfam�������ݿ��reads�İٷֱȣ�<br />
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

