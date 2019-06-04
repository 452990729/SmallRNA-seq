#!/usr/bin/perl

use strict;
use Bio::SeqIO;
die "Usage: perl $0 <fasta file> <prefix of output files> <number of seqs per file>.\nThis program splits user's fasta file to several files each of which contians given number seqs.\nSplitted files are named with user defined prefix. Bioperl module is needed.\n" unless @ARGV==3;

my $from = shift;
my $toprefix = shift;
my $seqs = shift;

my $in  = new Bio::SeqIO(-file  => $from);

my $count = 0;
my $fcount = 1;
#my $out = new Bio::SeqIO(-file => ">$toprefix\_$fcount.fasta", -format=>'fasta');
my $out;
while (my $seq = $in->next_seq) {
        if ($count % $seqs == 0) {
                $out = new Bio::SeqIO(-file => ">$toprefix\_$fcount.fasta", -format=>'fasta');
                $fcount++;
        }
        $out->write_seq($seq);
        $count++;
}
