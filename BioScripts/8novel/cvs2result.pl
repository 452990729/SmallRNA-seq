use warnings;
use strict;

use Bio::Perl;
use Bio::DB::Fasta;
use List::Util qw(reduce);
use FindBin '$Bin';
use Config::Tiny;


die "perl $0 cvs hp_db ma_db pos outdir" unless (@ARGV == 5);

my ($cvs_file,$hp_file,$ma_file,$pos_file,$outdir)=@ARGV;
if(-e "$hp_file.index"){
`rm -rf $hp_file.index`;
}
if(-e "$ma_file.index"){
`rm -rf $ma_file.index`;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $srnatools = $Config->{software}->{srnatoolscli};
my $hpin=Bio::DB::Fasta->new("$hp_file");
my $matin=Bio::DB::Fasta->new("$ma_file");

if(!(-e $outdir)){
	`mkdir -p $outdir`;
}
if(!(-e "$outdir/image")){
	`mkdir -p $outdir/image`;
}

`>$outdir/hairpin.fa`;
`>$outdir/hairpin_mature.fa`;
`>$outdir/mature.fa`;
`>$outdir/star.fa`;
`>$outdir/hairpin.pos`;

my %hash;
my %hash_ma;
my %hash_star;
my @hash_sample;
my %mat_sample;
my %star_sample;
open(TAG,$cvs_file);
open(CVS,">$outdir/miRNAs_expressed_novel.csv");
while(<TAG>){
	chomp;
	my @ary=split(/\t/,$_);
	if(!/^#/){
		if($ary[3]>0){
			open(HP,">$outdir/hairpin.fa_tmp");
			print HP ">$ary[2]\n";
			print HP $hpin->seq($ary[2])."\n";
			close(HP);
			open(MAT,">$outdir/mature.fa_tmp");
			print MAT ">$ary[0]\n";
			print MAT $matin->seq($ary[0])."\n";
			close(MAT);
			if(!defined($hash{$ary[2]})){
				`awk \'{if(\$1==\"$ary[2]\"){print}}\' $pos_file >>$outdir/hairpin.pos`;
				$hash{$ary[2]}=1;
				`cat $outdir/hairpin.fa_tmp >>$outdir/hairpin.fa`;
			}
			if(!defined($hash_ma{$ary[0]}) || !defined($hash_star{$ary[0]})){
				`cat $outdir/mature.fa_tmp >>$outdir/hairpin_mature.fa`;
			}
			if($ary[1]>0){
				print CVS join("\t",@ary),"\n";
				if($ary[0] !~ /\*$/){
					`$perlExec $srnatools/srna-tools.pl --tool hp_tool --longSeq $outdir/hairpin.fa_tmp --shortSeqs $outdir/mature.fa_tmp --out $outdir/image`;
					`mv $outdir/image/Structure_plot_bitmap.jpg $outdir/image/$ary[0]_$ary[2].jpg`;
					`mv $outdir/image/Structure_plot.pdf $outdir/image/$ary[0]_$ary[2].pdf`;
					if(!defined($hash_ma{$ary[0]})){
						my $marker=0;
						for my $i(0..$#hash_sample){
							if($ary[4+$i]>0){
								$mat_sample{$hash_sample[$i]}++;
								$marker=1;
							}
						}
						if($marker){$mat_sample{'Total'}++;}
						$hash_ma{$ary[0]}=1;
						`cat $outdir/mature.fa_tmp >>$outdir/mature.fa`;
					}
				}else{
					if(!defined($hash_star{$ary[0]})){
						my $marker=0;
						for my $i(0..$#hash_sample){
							if($ary[4+$i]>0){
								$star_sample{$hash_sample[$i]}++;
								$marker=1;
							}
						}
						if($marker){$star_sample{'Total'}++;}
						$hash_star{$ary[0]}=1;
						`cat $outdir/mature.fa_tmp >>$outdir/star.fa`;
					}
				}
			}
		}
	}else{
		my $sn=(scalar(@ary)-4)/2;
		for my $i(1..$sn){
			push @hash_sample,$ary[$i+3];
		}
		print CVS join("\t",@ary),"\n";
	}
}
close(CVS);

`rm -rf $outdir/hairpin.fa_tmp  $outdir/mature.fa_tmp  $outdir/image/legend.txt`;
open TOTAL,">$outdir/novel_miRNA.mapmat.stat";
my @join_sample;
my @join_sample_mat;
push @join_sample,"Total";
push @join_sample_mat,$mat_sample{"Total"};
for my $i(0..$#hash_sample){
	push @join_sample,$hash_sample[$i];
	push @join_sample_mat,$mat_sample{$hash_sample[$i]};
}
print TOTAL join("\t",@join_sample),"\n";
print TOTAL join("\t",@join_sample_mat),"\n";
close(TOTAL);

open TOTAL,">$outdir/novel_miRNA.mapstar.stat";
@join_sample=();
my @join_sample_star;
push @join_sample,"Total";
push @join_sample_star,$star_sample{"Total"};
for my $i(0..$#hash_sample){
        push @join_sample,$hash_sample[$i];
        push @join_sample_star,$star_sample{$hash_sample[$i]};
}
print TOTAL join("\t",@join_sample),"\n";
print TOTAL join("\t",@join_sample_star),"\n";
close(TOTAL);
