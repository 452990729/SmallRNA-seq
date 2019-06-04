use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;

die "perl $0 genome phase.predict phase.out\n" unless (@ARGV == 3);

my $genome = shift;
my $IN_dir=shift;
my $OUT_dir=shift;

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
$genome = get_ful_path($genome);

if(!(-e "$IN_dir")){
	die "directory : $IN_dir is not exist.\n";
}
if(!(-e "$IN_dir/locuslist.csv")|| !(-e "$IN_dir/srnas.txt")){
	die "locuslist.csv and srnas.txt not exist in $IN_dir, please make sure this is the output directory of srna-tools phasing analysis\n";
}

if(!(-e "$OUT_dir")){
	`mkdir -p $OUT_dir`;
}

open(IN1,"$IN_dir/locuslist.csv");
my %hash;
while(<IN1>){
	chomp;
	my @array=split/,/;
	if(/\d+/ && /^\S+\,\S+\,\S+\,\S+\,\S+\,/){
		if( $array[-1] < 1e-05 && $array[3] > 5 ){
			$hash{$array[0]}{'id'}++;
			$hash{$array[0]}{$hash{$array[0]}{'id'}}{'start'}=$array[1];
			$hash{$array[0]}{$hash{$array[0]}{'id'}}{'end'}=$array[2];
		}
	}
}
close(IN1);

open(OUT1,">$OUT_dir/novel_TAS.phased.fa");
my $in  = new Bio::SeqIO(-file => $genome);
while (my $seq = $in->next_seq()){
	my $start;
	my $end;
	if(defined($hash{$seq->id})){
		for my $i(1..$hash{$seq->id}{'id'}){
			$start=$hash{$seq->id}{$i}{'start'};
			$end=$hash{$seq->id}{$i}{'end'};
			if($start < 0){$start=1}
			if( $end > $seq->length ){$end=$seq->length}
			print OUT1 ">",$seq->id,":$start-$end\n",substr($seq->seq,$start-1,$end-$start+1),"\n";
		}
	}
}
close(OUT1);

sub get_ful_path{
        my $in=shift;
        if($in !~ /^\//){
                my $t=`pwd`;
                chomp($t);
                return "$t/$in";
        }else{
                return $in;
        }
}
