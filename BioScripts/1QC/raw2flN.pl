
#===============================================================================
=head1 Name

	pg_deadapter_se_v1.pl

=head1 Version

 Author:   Tian Shilin,tianshilin@novogene.cn
 Company:  NOVOGENE                                 
 Version:  1.0                                  
 Created:  09/15/2011 10:52:53 AM
 Modif: 06/26/2012 tianshilin@novogene.cn

=head1 Description
	
	-i|fq	(str)	input fq file(*.gz is OK)
	-l|len (f) sequence length
	-o|outdir	(str)	outdir   (default = ./)
	-nc|n-cutoff	(f)      in the PE reads (default = 0.1 )  the (N_num/len)>=$N is bad fq
	-q|qual	(int)    base quality type( default : 64 )
	-lq (int)   low quality(default : 5 )
	-qc|qual-cutoff	( default = 0.5 ) in the PE reads [ (num of quality <= $Qmin )/len ]>=$Qrate is bad fq
    -name Sample name( default: sample )

=head1 Usage
	
	perl pg_deadapter_se_v1.pl -i fq [ -l <N> -nc <N rate> -q <type quality> -lq <low quality> -qc <low quality rate> -o .]
	
=head1 Exmple

	perl pg_deadapter_se_v1.pl -i test.fq -nc 0.1 -q 33 -lq 5 -qc 0.5 -o ./ 

=cut

#===============================================================================

use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use List::Util qw(first max);
use FindBin '$Bin';
use Config::Tiny;

my $help;
my ($fq, $out_dir, $n_cutoff, $low_q, $low_qual_cutoff, $len, $type_q, $name);



GetOptions(
	"h|?|help"		=> \$help,
	"i|fq=s"		=> \$fq,
	"l|len=f"		=> \$len,
	"o|out_dir=s"		=> \$out_dir,
	"nc|n-cutoff=f"		=> \$n_cutoff,
	"q|qual=i"		=> \$type_q,
	"lq=i"			=> \$low_q,
	"qc|qual-cutoff=f"	=> \$low_qual_cutoff,
        "name=s"	        => \$name
);

die `pod2text $0` if $help;
die `pod2text $0` unless $fq;
my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $R_v2153 = $Config->{srnaenv}->{R_v2153};

#print STDERR "Filter begin at: ".`date`."\n";
$out_dir ||= ".";
$n_cutoff ||= 0.1;
$low_qual_cutoff ||= 0.5;
$type_q ||= 33;
$low_q ||= 5;
$name ||= 'sample';
my $q20 = $type_q + 19;
my $q30 = $type_q + 29;
my $low_qual = sprintf "\\%03o", ($low_q + $type_q);
##===================================================================#
$q20 = sprintf "\\%03o",$q20;
$q30 = sprintf "\\%03o",$q30;
my ($mark,$mark_ad);
if($fq =~ /\.gz$/){
	$mark = 1;
}else{
	$mark = 0;
}

#===================================================================#
# Read input fastq file and generate output.                        #
#===================================================================#
 
if($mark eq 1){
	open Fq, "<:gzip",$fq or die "ERROR: cann't ungzip $fq: $!";
}elsif($mark eq 0){
	open Fq,"$fq" || die $!;
}
#my $cleanfq = (split /\./,basename($fq))[0];
open Cleanfq,">","$out_dir/$name.clean.fq" or die $!;

my @reads_no = (0, 0);
my @N = (0, 0,);
my @GC = (0, 0);
my @low_qual = (0, 0);
my @q20 = (0, 0);
my @q30 = (0, 0);
my $discard_n = 0;
my $discard_low_qual = 0;
my @read_n = (0,0);
my (%base,%clean_base,%qual1,%qual2,%clean_qual1,%clean_qual2,%seq,%clean_seq);

while (1) {
	my (@q,$id);
	for(0..3){
		$q[$_] = <Fq>;
		if (defined $q[$_]){
			chomp ($q[$_]);
		}else{
			last;
		}
	}
	unless (defined $q[0]){
		last;
	}
	$reads_no[0] += length $q[1];
	$read_n[0] ++;
	my $N = $q[1] =~ tr/nN//;
	my $gc = $q[1] =~ tr/gcGC//;
	$N[0] += $N;
	$GC[0] += $gc;
	my $low_qual_num = eval "\$q[3] =~ tr/\0-$low_qual//";
	my $low_q20_num = eval "\$q[3] =~ tr/\0-$q20//";
	my $low_q30_num = eval "\$q[3] =~ tr/\0-$q30//";

	$low_qual[0] += $low_qual_num;
	$q20[0] += $low_q20_num;
	$q30[0] += $low_q30_num;
	$seq{"$q[1]"} = 1;

	my @base1 = split //,$q[1];
	my @qual1 = split //,$q[3];
	for (my $i=0;$i<(length $q[1]);$i++){
		$base{$i}{$base1[$i]} ++;
		$qual1{$i} += ((ord $qual1[$i]) - $type_q);
		$qual2{((ord $qual1[$i]) - $type_q)} += 1;
	}
	if (defined $n_cutoff && ((($n_cutoff*length $q[1]) <= $N))){
		$discard_n ++;
		next;
	}
	if (defined $low_qual_cutoff && ((($low_qual_cutoff*length $q[3]) <= $low_qual_num))){
		$discard_low_qual ++;
		next;
	}
#------------------------------------clean data stat---------------------------#
	for (my $i=0;$i<(length $q[1]);$i++){
		$clean_base{$i}{$base1[$i]} ++;
		$clean_qual1{$i} += ((ord $qual1[$i]) - $type_q);
		$clean_qual2{((ord $qual1[$i]) - $type_q)} += 1;
	}
#------------------------------------clean data stat---------------------------#

	$reads_no[1] += (length $q[1]);
	$read_n[1] ++;
	$N[1] += $N;
	$GC[1] += $gc;
	$low_qual[1] += $low_qual_num;
	$q20[1] += $low_q20_num;
	$q30[1] += $low_q30_num;
	$clean_seq{"$q[1]"} = 1;

	$q[0] =~ s/#\w+//g;
	my $read = join("\n",@q)."\n";
	print Cleanfq $read;
}
close Fq;
close Cleanfq;
my $key = $name;
#(my $key = $cleanfq) =~ s/\_1$//g;

$len ||= (int($reads_no[0]/($read_n[0]))+1);

open RGC,">$out_dir/raw\_$key.GC" || die $!;
open CGC,">$out_dir/clean\_$key.GC" || die $!;
open RQM,">$out_dir/raw\_$key.QM" || die $!;
open CQM,">$out_dir/clean\_$key.QM" || die $!;
open RQD,">$out_dir/raw_$key.QD" || die $!;
open CQD,">$out_dir/clean_$key.QD" || die $!;
my $error;
foreach my $base_site1(sort {$a<=>$b} keys %qual1){
	my $qm = ($qual1{$base_site1}/$read_n[0]);
	my $index = 0 -($qm/10);
	my $error_rate = (10**$index)*100;
	$error += $error_rate;
	print RQM $base_site1,"\t",$qm,"\t",$error_rate,"\n";
}
my $base_qual_max = max (keys %qual2);
foreach my $base_qual(sort {$a<=>$b} keys %qual2){
	print RQD $base_qual,"\t",$qual2{$base_qual},"\n";
}

foreach my $base_site2(sort {$a<=>$b} keys %base){
	print RGC $base_site2,"\t";
	foreach ("A","T","G","C","N"){
		$base{$base_site2}{$_} = 0 if (!$base{$base_site2}{$_});
		print RGC $_,"\t",$base{$base_site2}{$_},"\t",sprintf("%.2f",$base{$base_site2}{$_}*100/$read_n[0]),"\t"; 
	}
	print RGC "\n";
}

#------------------------------------clean data stat---------------------------#
my $clean_error;
foreach my $clean_base_site1(sort {$a<=>$b} keys %clean_qual1){
	my $clean_qm = ($clean_qual1{$clean_base_site1}/$read_n[1]);
	my $clean_index = 0 -($clean_qm/10);
	my $clean_error_rate = (10**$clean_index)*100;
	$clean_error += $clean_error_rate;
	print CQM $clean_base_site1,"\t",$clean_qm,"\t",$clean_error_rate,"\n";
}
my $clean_base_qual_max = max (keys %qual2);
foreach my $clean_base_qual(sort {$a<=>$b} keys %clean_qual2){
	print CQD $clean_base_qual,"\t",$clean_qual2{$clean_base_qual},"\n";
}

foreach my $clean_base_site2(sort {$a<=>$b} keys %clean_base){
	print CGC $clean_base_site2,"\t";
	foreach ("A","T","G","C","N"){
		$clean_base{$clean_base_site2}{$_} = 0 if (!$clean_base{$clean_base_site2}{$_});
		print CGC $_,"\t",$clean_base{$clean_base_site2}{$_},"\t",sprintf("%.2f",$clean_base{$clean_base_site2}{$_}*100/$read_n[1]),"\t"; 
	}
	print CGC "\n";
}

#------------------------------------clean data stat---------------------------#
#------------------------------------raw data plot ---------------------------#
my $C_GC =<< "END";
a<-read.table("$out_dir/raw\_$key.GC")
x<-a[,1]
y<-a[,4]
png("$out_dir/$name.raw_GC.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
plot(x,y,xlim=c(0,($len-1)),ylim=c(0,50),col="red",type="l",xlab="Position along reads",ylab="Percent",main="$name",lty=1,lwd=1.5)
p<-a[,7]
q<-a[,10]
s<-a[,13]
m<-a[,16]
lines(x,p,col="darkorange",type="l",lty=2,lwd=1.5)
lines(x,q,col="darkblue",type="l",lty=4,lwd=1.5)
lines(x,s,col="darkgreen",type="l",lty=5,lwd=1.5)
lines(x,m,col="cyan3",type="l",lty=6,lwd=1.5)
legend("topright",legend=c("A","T","G","C","N"),col=c("red","darkorange","darkblue","darkgreen","cyan3"),lty=c(1,2,4,5,6),bty='n')
dev.off()
END
my $C_QD =<< "END";
a<-read.table("$out_dir/raw\_$key.QD")
x<-a[,1]
y<-a[,2]
png("$out_dir/$name.raw_QD.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
plot(x,y,xlim=c(0,$base_qual_max),col="darkgreen",type="l",xlab="Quality score",ylab="Number of bases",main="$name")
axis(side=1,at=seq(from=0,to=$base_qual_max,by=5))
dev.off()
END
my $C_QM =<< "END";
a<-read.table("$out_dir/raw\_$key.QM")
x<-a[,1]
y<-a[,2]
z<-a[,3]
png("$out_dir/$name.raw_QM.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
plot(x,y,xaxt="n",xlim=c(0,($len-1)),ylim=c(0,40),col="darkgreen",type="p",pch=1,cex=1.5,xlab="Position along reads",ylab="Quality",main="$name")
axis(side=1,at=seq(from=0,to=($len-1),by=20))
abline(h=20,col="darkblue",lty=3)
abline(v=20,col="darkblue",lty=3)
abline(v=40,col="darkblue",lty=3)
abline(v=60,col="darkblue",lty=3)
abline(v=80,col="darkblue",lty=3)
abline(v=100,col="darkblue",lty=3)
abline(v=120,col="darkblue",lty=3)
abline(v=140,col="darkblue",lty=3)
abline(v=160,col="darkblue",lty=3)
abline(v=180,col="darkblue",lty=3)
png("$out_dir/$name.raw_Error.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
plot(x,z,xaxt="n",xlim=c(0,($len-1)),col="darkgreen",type="h",xlab="Position along reads",ylab="% Error-rate",main="$name")
axis(side=1,at=seq(from=0,to=($len-1),by=20))
abline(v=20,col="darkblue",lty=3)
abline(v=40,col="darkblue",lty=3)
abline(v=60,col="darkblue",lty=3)
abline(v=80,col="darkblue",lty=3)
abline(v=100,col="darkblue",lty=3)
abline(v=120,col="darkblue",lty=3)
abline(v=140,col="darkblue",lty=3)
abline(v=160,col="darkblue",lty=3)
abline(v=180,col="darkblue",lty=3)
dev.off()
END
#------------------------------------raw data plot ---------------------------#
open R,"|$R_v2153/R --vanilla --slave" or die $!;
print R $C_GC;
print R $C_QM;
print R $C_QD;
close R;

#------------------------------------clean data plot ---------------------------#
$C_GC =<< "END";
a<-read.table("$out_dir/clean\_$key.GC")
x<-a[,1]
y<-a[,4]
png("$out_dir/$name.clean_GC.png",type="cairo-png",,width=480*4,height=480*4,res=72*4)  
plot(x,y,xlim=c(0,($len-1)),ylim=c(0,50),col="red",type="l",xlab="Position along reads",ylab="Percent",main="$name",lty=1,lwd=1.5)
p<-a[,7]
q<-a[,10]
s<-a[,13]
m<-a[,16]
lines(x,p,col="darkorange",type="l",lty=2,lwd=1.5)
lines(x,q,col="darkblue",type="l",lty=4,lwd=1.5)
lines(x,s,col="darkgreen",type="l",lty=5,lwd=1.5)
lines(x,m,col="cyan3",type="l",lty=6,lwd=1.5)
legend("topright",legend=c("A","T","G","C","N"),col=c("red","darkorange","darkblue","darkgreen","cyan3"),lty=c(1,2,4,5,6),bty='n')
dev.off()
END
$C_QD =<< "END";
a<-read.table("$out_dir/clean\_$key.QD")
x<-a[,1]
y<-a[,2]
png("$out_dir/$name.clean_QD.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
plot(x,y,xlim=c(0,$clean_base_qual_max),col="darkgreen",type="l",xlab="Quality score",ylab="Number of bases",main="$name")
axis(side=1,at=seq(from=0,to=$base_qual_max,by=5))
dev.off()
END
$C_QM =<< "END";
a<-read.table("$out_dir/clean\_$key.QM")
x<-a[,1]
y<-a[,2]
z<-a[,3]
png("$out_dir/$name.clean_QM.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
plot(x,y,xaxt="n",xlim=c(0,($len-1)),ylim=c(0,40),col="darkgreen",type="p",pch=1,cex=1.5,xlab="Position along reads",ylab="Quality",main="$name")
axis(side=1,at=seq(from=0,to=($len-1),by=20))
abline(h=20,col="darkblue",lty=3)
abline(v=20,col="darkblue",lty=3)
abline(v=40,col="darkblue",lty=3)
abline(v=60,col="darkblue",lty=3)
abline(v=80,col="darkblue",lty=3)
abline(v=100,col="darkblue",lty=3)
abline(v=120,col="darkblue",lty=3)
abline(v=140,col="darkblue",lty=3)
abline(v=160,col="darkblue",lty=3)
abline(v=180,col="darkblue",lty=3)
png("$out_dir/$name.clean_Error.png",type="cairo-png",width=480*4,height=480*4,res=72*4)
plot(x,z,xaxt="n",xlim=c(0,($len-1)),col="darkgreen",type="h",xlab="Position along reads",ylab="% Error-rate",main="$name")
axis(side=1,at=seq(from=0,to=($len-1),by=20))
abline(v=20,col="darkblue",lty=3)
abline(v=40,col="darkblue",lty=3)
abline(v=60,col="darkblue",lty=3)
abline(v=80,col="darkblue",lty=3)
abline(v=100,col="darkblue",lty=3)
abline(v=120,col="darkblue",lty=3)
abline(v=140,col="darkblue",lty=3)
abline(v=160,col="darkblue",lty=3)
abline(v=180,col="darkblue",lty=3)
dev.off()
END
#------------------------------------clean data plot ---------------------------#

open R,"| $R_v2153/R --vanilla --slave" or die $!;
print R $C_GC;
print R $C_QM;
print R $C_QD;
close R;


open Stat,">$out_dir/$key.stat" || die $!;
open RawStat,">$out_dir/${key}_raw.stat" || die $!;
my $rate = sprintf("%.2f%%",$reads_no[1]*100/$reads_no[0]);
my @low_qual_rate = (0,0);
$low_qual_rate[0] = sprintf("%.2f%%",$low_qual[0]*100/$reads_no[0]);
$low_qual_rate[1] = sprintf("%.2f%%",$low_qual[1]*100/$reads_no[1]);
my @q20_rate = (0,0);
$q20_rate[0] = sprintf("%.2f%%",(100-$q20[0]*100/$reads_no[0]));
$q20_rate[1] = sprintf("%.2f%%",(100-$q20[1]*100/$reads_no[1]));
my @q30_rate = (0,0);
$q30_rate[0] = sprintf("%.2f%%",(100-$q30[0]*100/$reads_no[0]));
$q30_rate[1] = sprintf("%.2f%%",(100-$q30[1]*100/$reads_no[1]));
my @GC_rate = (0,0);
$GC_rate[0] = sprintf("%.2f%%",$GC[0]*100/($reads_no[0]-$N[0]));
$GC_rate[1] = sprintf("%.2f%%",$GC[1]*100/($reads_no[1]-$N[1]));
my @N_rate = (0,0);
$N_rate[0] = sprintf("%.2f%%",$N[0]*100/$reads_no[0]);
$N_rate[1] = sprintf("%.2f%%",$N[1]*100/$reads_no[1]);
my @dup_rate = (0,0);
my $uniq = scalar (keys %seq);
$dup_rate[0]=sprintf("%.2f%%",($read_n[0]-$uniq)*100/$read_n[0]);
$uniq = scalar (keys %clean_seq);
$dup_rate[1]=sprintf("%.2f%%",($read_n[1]-$uniq)*100/$read_n[1]);
my $DataSize = (sprintf "%.3f",$reads_no[0] * 1e-9)."G";
my $ErrorRate = sprintf("%.2f%%",$error/100);

print Stat "Type\tRaw data\tClean data\n";
print Stat "Number of Reads:\t$read_n[0]\t$read_n[1]\n";
print Stat "Data Size:\t$reads_no[0]\t$reads_no[1]\($rate\)\n";
print Stat "N of fq:\t$N_rate[0]\t$N_rate[1]\n";
print Stat "Low qual base(q<=5) of fq:\t$low_qual_rate[0]\t$low_qual_rate[1]\n";
print Stat "Q20 of fq:\t$q20_rate[0]\t$q20_rate[1]\n";
print Stat "Q30 of fq:\t$q30_rate[0]\t$q30_rate[1]\n";
print Stat "GC of fq:\t$GC_rate[0]\t$GC_rate[1]\n";
print Stat "Error of fq:\t",sprintf("%.2f%%",$error/100),"\t",sprintf("%.2f%%",$clean_error/100),"\n";
print Stat "Duplication Rate of fq:\t$dup_rate[0]\t$dup_rate[1]\n";
print Stat "Discard Reads related to N:\t$discard_n\n";
print Stat "Discard Reads related to low qual:\t$discard_low_qual\n";
print Stat "Reads Classification:\t$read_n[0]\t$read_n[1]\t$discard_n\t$discard_low_qual\n";

print RawStat "Sample\tReads\tBases\tError rate\tQ20\tQ30\tGC content\n";
print RawStat "$key\t$read_n[0]\t$DataSize\t$ErrorRate\t$q20_rate[0]\t$q30_rate[0]\t$GC_rate[0]";
system("cut -f 1,3 $out_dir/raw_${key}.QM > $out_dir/${key}.raw_Error.txt");

#print STDERR "Filter finaly at:\t".`date`."\n";

#Run cut the first several base
system("awk \'{if(NR%2==1){print}else{print substr(\$1,0+1,50)}}\' $out_dir\/${key}.clean.fq > $out_dir\/${key}.hq.cut.fq && \\
      rm -rf $out_dir\/${key}.clean.fq");
