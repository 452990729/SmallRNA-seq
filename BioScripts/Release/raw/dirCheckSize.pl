#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
my ($s,$c) = ("checkSize.xls");
GetOptions("c"=>\$c,"s:s"=>\$s);
@ARGV || die"
Usage1(make checkSize file): perl CheckSize.pl <dir/> 
    -s <file>       set output checkSize file name, defualt=checkSize.xls
Usage2(check dir completeness): perl CheckSize.pl -c <dirs/ | dir.list>
    -s <file>       set input checkSize file name in dir/ or each dir at dir.list, defualt=checkSize.xls
    dir.list form:  dir/\tcheckSize_file_name in each line, checkSize file name default set by -s\n\n";
if($c){
    my @dirs;
    for my $d(@ARGV){
        if(-f $d && -s $d){
            open IN,$d || die$!;
            while(<IN>){
                /\S/ || next;
                my @l = split;
                push @dirs,[$l[0],$l[1]||$s];
            }
            close IN;
        }else{
            push @dirs,[$d,$s];
        }
    }
    for my $p(@dirs){
        my $sign = check(@$p);
        print "$p->[0]\t$sign\n";
    }
}else{
    my @files;
	my $f = shift;
	(-d $f) || die"Note: input ARGV should be directory!\n";
	chdir $f;    
	my @lnfs;
    for (`find ./ -type l`){chomp;push @lnfs,$_;}
    pureln(\@files,\@lnfs);
	for (`find ./ -type f`){
	    chomp;
	    ($_ eq "./$s") && next;
	    push @files,$_;
	}
	open OUT,">$s" || die$!;
	for my $f(@files){
	    my $af = (-l $f) ? abs_path($f) : $f;
	    my $size = (-s $f) ? (split/\s+/,`ls -l $af`)[4] : 0;
	    ($size == 0) && (-l $f) && die"erro link: $f to $af\n";
	    print OUT "$size\t$f\n";
	}
	close OUT;
}
sub pureln{
    my ($files,$lnfs) = @_;
    @$lnfs || return(0);
    my @newln;
    for my $f(@$lnfs){
        my $af = (-l $f) ? abs_path($f) : $f;
        if(-d $af){
            chomp(my @f1 = `find $f/ -type l`);
            chomp(my @f2 = `find $f/ -type f`);
            @f1 && pureln($files,\@f1);
            @f2 && (push @$files,@f2);
        }else{
            push @$files,$f;
        }
    }
}
sub check{
    my ($dir,$sizef) = @_;
    (-d $dir) || return("can't find directory");
    (-s "$dir/$sizef") || return("can't find check file: $sizef");
    open IN,"$dir/$sizef" || die$!;
    while(<IN>){
        my @l = split;
        my $size = (-s "$dir/$l[1]") ? (split/\s+/,`ls -l $dir/$l[1]`)[4] : 0;
        ($size == $l[0]) || return("error: $l[1]");
    }
    close IN;
    return "OK";
}
