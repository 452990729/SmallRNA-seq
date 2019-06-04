#!/bin/perl -w
use strict;
use Getopt::Long;
use Config::Tiny;

my ($help, $query, $subject, $evalue, $maxTarget, $blastout, $gff, $transout);
GetOptions(
	"h|help"	=>\$help,
	"q|query=s"	=>\$query,
	"s|subject=s" =>\$subject,
	"e|evalue:f"	=>\$evalue,
	"m|maxTarget:i"	=>\$maxTarget,
	"bo|blastout=s"	=>\$blastout,
	"g|gff=s"		=>\$gff,
	"to|transout=s"	=>\$transout
);

my $usage=<<END;
-------------------------------------------------------------------------
	perl $0 -q query.fa -s subject.fa -e 1e-5 -m 3 -bo blastout -gff *.gff  -to transout

	-h|help		help
	-q|query	query for blastn
	-s|subject	subject for blastn
	-e|evalue	evalue for blastn
	-m|maxTarget	max_target_seqs for blastn
	-bo|blastout	blastn outputfile
	-g|gff		gff file
	-to|transout	tran NAT outputfile
-------------------------------------------------------------------------
END
die $usage if ($help or !$query or !$subject);

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $blast = $Config->{software}->{blast2};

## ------------------------------------- Blast -------------------------------------------------
$evalue||=1e-5;
$maxTarget||=3;
if( $blastout)
{
	`$blast/makeblastdb -dbtype nucl -in $subject`; 
	`$blast/blastn -db $subject -query $query -evalue $evalue -max_target_seqs $maxTarget -num_threads 4 -strand minus -out $blastout -outfmt "6 qseqid qlen sseqid slen length qstart qend sstart send evalue bitscore pident mismatch gaps"`;
} 

## ------------------------------------- Subroutine -------------------------------------------------
my ($Q_chr, $Q_strand, $Q_chrStart, $Q_chrEnd, $Q_id, $S_chr, $S_strand, $S_chrStart, $S_chrEnd, $S_id);
sub prefix()
{
	print O "$Q_chr\t$Q_strand\t$Q_chrStart\t$Q_chrEnd\t$Q_id\t$S_chr\t$S_strand\t$S_chrStart\t$S_chrEnd\t$S_id\t";
}

my ($trans_len, $Q_trans_chrStart, $Q_trans_chrEnd, $S_trans_chrStart, $S_trans_chrEnd, $Q_trans_geneStart, $Q_trans_geneEnd, $S_trans_geneStart, $S_trans_geneEnd);
sub suffix()
{
	print O "\t$trans_len\t$Q_trans_chrStart\t$Q_trans_chrEnd\t$S_trans_chrStart\t$S_trans_chrEnd\t$Q_trans_geneStart\t$Q_trans_geneEnd\t$S_trans_geneStart\t$S_trans_geneEnd\n";
}

## ------------------------------------- Main body -------------------------------------------------
open G, $gff or die "Can't open $gff: $!";  # gff file used to generate the chr postion of query/subject gene 
my %hash;
while(<G>)
{
    chomp;
    if(/\S+\t\S+\tgene\t/){
         my ($chr, $source, $type, $start, $end, $score, $strand, $frame, $attr)=split;
         my $id=$1 if ($attr=~/ID=([^;]+);/);
        $hash{$id}=[($chr, $start, $end, $strand)];
    }   
}
close G;

open O, ">$transout" or die "Can't output $transout: $!";
print O "Q_chr\tQ_strand\tQ_chrStart\tQ_chrEnd\tQ_id\tS_chr\tS_strand\tS_chrStart\tS_chrEnd\tS_id\tTrans\tTrans_len\tQ_trans_chrStart\tQ_trans_chrEnd\tS_trans_chrStart\tS_trans_chrEnd\tQ_trans_geneStart\tQ_trans_geneEnd\tS_trans_geneStart\tS_trans_geneEnd\n";

open B, "<$blastout" or die "Can't open $blastout: $!";   # Parse the blastn result to generate the NAT-trans pairs!
my %hash2;
while(<B>)
{
    chomp;
	#AT1G01010       2268    ATMG01400       30949358        99      1913    2010    2870    2776    3e-23    111    87.88   7       5
	my($qid, $qlen, $sid, $slen, $Alen, $qstart, $qend, $send, $sstart, $evalue, $score, $pidentiy, $mismatch, $gaps)=split; #Note: send->sstart
	($trans_len, $Q_trans_geneStart, $Q_trans_geneEnd, $S_trans_geneStart, $S_trans_geneEnd)=($Alen, $qstart, $qend, $sstart, $send);
	if(!defined $hash2{$qid} || !defined $hash2{$sid})
	{
		$hash2{$qid}=1;	$hash2{$sid}=1;
	    if($gaps!=0 && $gaps/$Alen<0.1 && ($Alen>=0.5*$qlen || $Alen>=0.5*$slen))
    	{
    		#if(defined $hash{$Q_id} && defined $hash{$S_id})
    		if(defined $hash{$qid} && defined $hash{$sid})
    		{   
				($Q_chr, $Q_chrStart, $Q_chrEnd, $Q_strand, $Q_id)=(@{$hash{$qid}}, $qid);
				($S_chr, $S_chrStart, $S_chrEnd, $S_strand, $S_id)=(@{$hash{$sid}}, $sid);

				$Q_trans_chrStart=$qstart+$Q_chrStart-1;
				$Q_trans_chrEnd=$qend+$Q_chrStart-1;  #Note: qend+Q_chrStart-1

				$S_trans_chrStart=$sstart+$S_chrStart-1;
				$S_trans_chrEnd=$send+$S_chrStart-1;  # end+start

				&prefix; print O "HC"; &suffix;
    		}   
    	}
    	elsif($gaps==0 && $trans_len>=100)
    	{
    		if(defined $hash{$qid} && defined $hash{$sid})
    		{   
				($Q_chr, $Q_chrStart, $Q_chrEnd, $Q_strand, $Q_id)=(@{$hash{$qid}}, $qid);
				($S_chr, $S_chrStart, $S_chrEnd, $S_strand, $S_id)=(@{$hash{$sid}}, $sid);

				$Q_trans_chrStart=$qstart+$Q_chrStart-1;
				$Q_trans_chrEnd=$qend+$Q_chrStart-1;  #Note: qend+Q_chrStart-1

				$S_trans_chrStart=$sstart+$S_chrStart-1;
				$S_trans_chrEnd=$send+$S_chrStart-1;  # end+start

				&prefix; print O "100nt"; &suffix;
    		}   
    	}
	}
}
close B;

