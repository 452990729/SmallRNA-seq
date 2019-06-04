use warnings;
use strict;
use Bio::SeqIO;

die "Usage: perl $0 id.list db.fasta out.fasta\nfetch the list of id from db.fasta to out.fasta\n" unless (scalar(@ARGV)==3);

my %id;
open(IN,"$ARGV[0]");
while(<IN>){
	chomp;
	/(\S+)/;
	$id{$1}=1;
#	print $1,"\n";
}

my $in  = new Bio::SeqIO(-file => $ARGV[1]);
my $out = new Bio::SeqIO(-file => ">$ARGV[2]", -format=>'fasta');
while (my $seq = $in->next_seq()){
	if(defined($id{$seq->id})){
		$out->write_seq($seq);
	}
}
