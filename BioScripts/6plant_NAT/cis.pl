#!/bin/perl -w
use strict;
use Getopt::Long;
my ($help, $gff, $outdir, $outfile);
GetOptions(
	"h|help"		=>\$help,
	"g|gff=s"		=>\$gff,
	"od|outdir=s"	=>\$outdir,
	"o|outfile=s"	=>\$outfile
);
my $usage=<<END;
---------------------------------------------------------
	perl $0 -g *.gff -od outdir -o outfile 
	-h|help		help
	-g|gff		gff file
	-od|outdir	outputdir
	-o|outfile	outputfile
---------------------------------------------------------
END
die $usage if ($help or !$gff);

##---------------------------------------- GFF Preprocess -----------------------------------------------------------------------
`awk 'sub(/;.+/, "") && \$3=="gene" && \$7=="+" && OFS="\t"{print \$1,\$4,\$5,\$7,\$9}' $gff|sed -e 's/ID=//g'|sort -k1,1 -k2,2n >$outdir/gff_pos.sort`;
`awk 'sub(/;.+/, "") && \$3=="gene" && \$7=="-" && OFS="\t"{print \$1,\$4,\$5,\$7,\$9}' $gff|sed -e 's/ID=//g'|sort -k1,1 -k2,2n >$outdir/gff_neg.sort`;

##---------------------------------------- Subroutine  -----------------------------------------------------------------------
my ($p, $n);
my ($pchr, $phead, $ptail, $pstrand, $pid, $nchr, $nhead, $ntail, $nstrand, $nid);
sub upP()
{
	$p=<P>;
	last unless defined $p;
	chomp $p;
	($pchr, $phead, $ptail, $pstrand, $pid)=split/\t/,$p;
}

sub upN()
{
	$n=<N>;
	last unless defined $n;
	chomp $n;
	($nchr, $nhead, $ntail, $nstrand, $nid)=split/\t/,$n;
}

sub prefix()
{
	print O "$pchr\t+\t$phead\t$ptail\t$pid\t$nchr\t-\t$nhead\t$ntail\t$nid\t";
}

my ($cis_len, $P_cis_chrStart, $P_cis_chrEnd, $N_cis_chrStart, $N_cis_chrEnd, $P_cis_geneStart, $P_cis_geneEnd, $N_cis_geneStart, $N_cis_geneEnd);
sub suffix()
{
	print O "\t$cis_len\t$P_cis_chrStart\t$P_cis_chrEnd\t$N_cis_chrStart\t$N_cis_chrEnd\t$P_cis_geneStart\t$P_cis_geneEnd\t$N_cis_geneStart\t$N_cis_geneEnd\n";
}

##---------------------------------------- Main body  -----------------------------------------------------------------------
open P, "<$outdir/gff_pos.sort" or die "Can't open: $!";
chomp ($p=<P>);
($pchr, $phead, $ptail, $pstrand, $pid)=split/\t/,$p;

open N, "<$outdir/gff_neg.sort"  or die "Can't open: $!";
chomp ($n=<N>);
($nchr, $nhead, $ntail, $nstrand, $nid)=split/\t/,$n;

open O, ">$outdir/$outfile";
print O "P_chr\tP_strand\tP_chrStart\tP_chrEnd\tP_id\tN_chr\tN_strand\tN_chrStart\tN_chrEnd\tN_id\tCis\tCis_len\tP_cis_chrStart\tP_cis_chrEnd\tN_cis_chrStart\tN_cis_chrEnd\tP_cis_geneStart\tP_cis_geneEnd\tN_cis_geneStart\tN_cis_geneEnd\n";

while (1)
{
	if($pchr lt $nchr)
	{
		&upP;
	}
	elsif ($pchr gt $nchr)
	{
		&upN;
	}
	else
	{
		if ($phead-$ntail>100)
		{
			&upN;
		}
		elsif ($phead-$ntail<=100 && $phead-$ntail>0)
		{
			$cis_len=$phead-$ntail+1; 
			($P_cis_chrStart, $P_cis_chrEnd)=($ntail, $phead); 
			($N_cis_chrStart, $N_cis_chrEnd)=($ntail, $phead); 
			($P_cis_geneStart, $P_cis_geneEnd)=(-1,-1);
			($N_cis_geneStart, $N_cis_geneEnd)=(-1,-1);
			&prefix; print O "Nearby head to head"; &suffix; 
			&upN;
		}
		elsif ($nhead<$phead && $ntail>=$phead && $ntail<$ptail)  # Note: ntail=phead !!
		{
			$cis_len=$ntail-$phead+1;
			($P_cis_chrStart, $P_cis_chrEnd)=($phead, $ntail);
			($N_cis_chrStart, $N_cis_chrEnd)=($phead, $ntail);
			($P_cis_geneStart, $P_cis_geneEnd)=(1, $cis_len);
			($N_cis_geneStart, $N_cis_geneEnd)=(1, $cis_len);
			&prefix; print O "Head to head"; &suffix;
			&upN;
		}
		elsif ($nhead<=$phead && $ntail>=$ptail)
		{
			$cis_len=$ptail-$phead+1;
			($P_cis_chrStart, $P_cis_chrEnd)=($phead, $ptail);
			($N_cis_chrStart, $N_cis_chrEnd)=($phead, $ptail);
			($P_cis_geneStart, $P_cis_geneEnd)=(1, $cis_len);
			my $negG_len=$ntail-$nhead+1;
			($N_cis_geneStart, $N_cis_geneEnd)=(($negG_len-($ptail-$nhead)), ($negG_len-($ptail-$nhead)+$cis_len));
			&prefix; print O "Containing"; &suffix;
			&upP;
		}
		elsif ($nhead >=$phead && $ntail<=$ptail)
		{
			$cis_len=$ntail-$nhead+1;
			($P_cis_chrStart, $P_cis_chrEnd)=($nhead, $ntail);
			($N_cis_chrStart, $N_cis_chrEnd)=($nhead, $ntail);
			my $posG_len=$ptail-$phead;
			($P_cis_geneStart, $P_cis_geneEnd)=(($posG_len-($ptail-$nhead)), ($posG_len-($ptail-$nhead)+$cis_len+1));
			($N_cis_geneStart, $N_cis_geneEnd)=(1, $cis_len);
			&prefix; print O "Containing"; &suffix;
			&upN;
		}
		elsif ($nhead>$phead && $nhead<=$ptail && $ntail>$ptail) #Note: nhead=ptail !!
		{
			$cis_len=$ptail-$nhead+1;
			($P_cis_chrStart, $P_cis_chrEnd)=($nhead, $ptail);	
			($N_cis_chrStart, $N_cis_chrEnd)=($nhead, $ptail);	
			my $posG_len=$ptail-$phead; 
			my $negG_len=$ntail-$nhead;
			($P_cis_geneStart, $P_cis_geneEnd)=(($posG_len-$cis_len), ($posG_len+1));
			($N_cis_geneStart, $N_cis_geneEnd)=(($negG_len-$cis_len), ($negG_len+1));
			&prefix; print O "Tail to tail"; &suffix;
			&upP;
		}
		elsif ($nhead-$ptail<=100 && $nhead-$ptail>0)
		{
			$cis_len=$nhead-$ptail+1;
			($P_cis_chrStart, $P_cis_chrEnd)=($ptail, $nhead);
			($N_cis_chrStart, $N_cis_chrEnd)=($ptail, $nhead);
			($P_cis_geneStart, $P_cis_geneEnd)=(-1, -1);
			($N_cis_geneStart, $N_cis_geneEnd)=(-1, -1);
			&prefix; print O "Nearby tail to tail"; &suffix;
			&upP;
		}
		elsif ($nhead-$ptail>100) 
        {   
           &upP;
        }   
	}
}
close P; close N; close O;
