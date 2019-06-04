use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;
use FindBin '$Bin';
use Config::Tiny;

die "perl $0 genome known_TAS.fa(TAIR10_rice.TAS.fna)" unless (@ARGV == 2);

my $genome = shift;
my $TAS_db=shift;

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $blast = $Config->{software}->{blast2};
my $pwd = `pwd`;
$pwd =~ s/\s+//g;
$genome = get_ful_path($genome);
$TAS_db = get_ful_path($TAS_db);

if(!(-e "known_TAS")){
	`mkdir -p known_TAS`;
}

if(!(-e "$genome.nhd")){
	if(!(-e "known_TAS/ref")){
		`mkdir -p known_TAS/ref`;
	}
	my $file=basename($genome);
	if(!(-e "known_TAS/ref/$file")){
		`ln -sf $genome known_TAS/ref/$file`;
	}
	if(!(-e "known_TAS/ref/$file.nhd")){
		`$blast/makeblastdb -in known_TAS/ref/$file -dbtype nucl -parse_seqids -hash_index`;
	}
	$genome="known_TAS/ref/$file";
}
if(!(-e "known_TAS/tmp")){
	`mkdir -p known_TAS/tmp`;
}
`$blast/blastn -query $TAS_db -db $genome -evalue 0.01 -outfmt 6 -out known_TAS/tmp/known_TAS.blastn`;
=head
open(FILE,"known_TAS/tmp/known_TAS.blastn");
`>known_TAS/known_TAS.blastn.fa`;
while(<FILE>){
        chomp;
        my @array=split/\t/;
	my $id = $array[1];
	my $length=`$blast_bin_dir/blastdbcmd -db $genome -dbtype nucl -entry $id |awk '{if(NR>1){len+=length(\$0)}}END{print len}'`;
	chomp($length);
	if($array[8]<$array[9]){
		my ($start,$end)=($array[8]-300,$array[9]+300);
		if($start<0){$start=1}
		if($end>$length*1){$end=$length}
		`$blast_bin_dir/blastdbcmd -db $genome -dbtype nucl -entry $id -range $start-$end|sed \'s/lcl|//g\' >>known_TAS/known_TAS.blastn.fa`;
	}else{
		my ($start,$end)=($array[9]-300,$array[8]+300);
		if($start<0){$start=1}
		if($end>$length*1){$end=$length}
		`$blast_bin_dir/blastdbcmd -db $genome -dbtype nucl -entry $id -range $start-$end -strand minus|sed \'s/lcl|//g\' >>known_TAS/known_TAS.blastn.fa`;
	}
}
=cut
my($id, %hash);
open G, $genome or die $!;
open O, ">known_TAS/known_TAS.blastn.fa" or die $!;
while(<G>)
{
	chomp;
	if(/>(\S+)/)
	{
		$id=$1;
	}
	else
	{
		$hash{$id}.=$_;
	}
}
close G;

open B, "known_TAS/tmp/known_TAS.blastn" or die $!;
while(<B>)
{
	chomp;
#TAS1b	1	100.00	839	0	0	1	839	18550042	18549204	0.0	1550
#TAS1b	2	87.22	407	33	9	1	396	11722466	11722068	5e-124	 446
#TAS1b	2	82.32	413	50	14	1	399	16538277	16537874	3e-91	 337
	my ($qid, $sid, $identity, $len, $mis, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $score)=split/\t/,$_;
	if(defined $hash{$sid})
	{
		if($sstart<$send)
		{
			my $seq=substr($hash{$sid}, $sstart-300-1, $send+300-$sstart+1);
			my $low=$sstart-300; my $up=$send+300;
			print O ">$sid:$low-$up\n$seq\n";
		}
		else
		{
			my $seq=substr($hash{$sid}, $send-300-1, $sstart+300-$send+1);
			$seq=~tr/ATCGatcg/TAGCtagc/;
			$seq=reverse $seq;
			my $up=$send-300; my $low=$sstart+300; # Note: low>up, mean minus strand 
			print O ">$sid:$low-$up\n$seq\n";
		}
	}
}	
close B;
close O;

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

sub min_max{
	my ($x,$y)=@_;
	if($x>$y){
		my $tmp=$x;
		$x=$y;
		$y=$x;
	}
	return ($x,$y);
}
