use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($help,$gtf,$tp,$out);
GetOptions(
	"h|?|help" => \$help,
	"gtf=s" => \$gtf,
	"tp=s" => \$tp,
	"o=s" => \$out,
);

open GTF,"<$gtf";
open TP,"<$tp";
open OUT,">$out";

my %hash;
my $gene_id;
my $transcript_id;
print OUT "miRNA\ttarget_mRNA\ttarget_gene\n";
while (<GTF>){
	chomp;
	if (/gene_id "(\S+)";/){
		$gene_id=$1;
	}
	if (/transcript_id "(\S+)";/){
		$transcript_id=$1;
	}
		$hash{$transcript_id}=$gene_id;
	}

while (<TP>){
	chomp;
        next if $_ =~ m/^#/;
	my @tmp=split /\t/;
	if (exists $hash{$tmp[1]}){
		print OUT "$tmp[0]\t$tmp[1]\t$hash{$tmp[1]}\n";
	}
}

