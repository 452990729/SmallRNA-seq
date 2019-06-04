#!/usr/bin/perl -w
use lib "/usr/local/Bioperl/lib/perl5/site_perl/5.8.3";
use Bio::DB::Fasta;
use strict;
my $seq_file = shift;
chomp $seq_file;
my $seqinx = Bio::DB::Fasta->new("$seq_file");
exit 0;
