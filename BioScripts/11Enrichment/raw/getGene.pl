#!/PROJ/RNA/share/software/perl/bin/perl -w
use strict;
use Bio::Perl;
use Bio::DB::Fasta;
use Data::Dumper;

if (@ARGV < 3) {
	print "Usage: $0 ReferenceDir Transcript.gtf Output.fa\n";
	exit;
}
my ($referenceDir,$gtfile,$output) = @ARGV;
my (%gene,%trans);
my $db = Bio::DB::Fasta->new("$referenceDir");
open GTF,"<$gtfile" or die $!;
while (<GTF>) {
	chomp;
	my @line = split /\t/,$_;
	next if ($line[2] ne "exon");
	#print "$line[8]\n";
	next unless ($line[8] =~ /gene_id "(\S+)"; transcript_id "(\S+)"/);
	my $geneID = $1;
	my $transID = $2;
	#print "$geneID\t$transID\n";
	push(@{$gene{$geneID}},$transID);
	my $exonSeq= ($line[6] eq "+") ? $db->subseq($line[0],$line[3],$line[4]) : $db->subseq($line[0],$line[4],$line[3]);
        #print "$line[0]	$line[8]\n" if(!$exonSeq);
	push (@{$trans{$transID}},[$line[3],$exonSeq]);
}
close GTF;
foreach my $transID (keys %trans) {
	$trans{$transID} = &link($trans{$transID});
}
#print Dumper(\%trans);
#print Dumper(\%gene);
foreach my $geneID (keys %gene) {
	$gene{$geneID} = &getGeneSeq($gene{$geneID});
}
open OUT,">$output" or die $!;
foreach my $geneID (sort keys %gene) {
	print OUT ">$geneID\n$gene{$geneID}\n";
}
close OUT;
sub getGeneSeq {
	my @ids = @{$_[0]};
	my $maxLen = 0;
	my $string;
	foreach my $id (@ids) {
		my $len = length($trans{$id});
		if ($maxLen < $len) {
			$maxLen = $len;
			$string = $trans{$id};
		}
	}
	return $string;
}




sub link {
	my @array = @{$_[0]};
	@array = sort {$a->[0] <=> $b->[0]}@array;
	my $seq;
	foreach (@array) {
		$seq .= $_->[1];
	}
	return $seq;
}
