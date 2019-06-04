use strict;
use warnings;


open IN,"<$ARGV[0]";
open DEL,">$ARGV[1]";

my @tmp;

while(<IN>){
	chomp;
	@tmp=split / /;
	if($tmp[-1] =~/tar.gz|md5sum.txt$|.R$|.py$|.pl$|.sh$|mapfile|shell/){
		next;
	}else{
		print DEL "rm $tmp[-1]\n";
	}
}
my $pwd=`pwd`;
my @temp;
@temp=split / /,$pwd;
my @tmpp;
@tmpp=split /\//,$temp[-1];
my $pro=$tmpp[-2];

#print DEL "perl /PUBLIC/source/RNA/smallRNA/version3/dir_process/change_step3_4.pl /ifs/TJPROJ3/RNA/WORK/project_status/$pro*";

