use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Bio::SeqIO;
use Config::Tiny;


sub usage{
print STDERR <<USAGE;
=========================================================================
Usage: perl $0 [options]
example: echo "perl $0 -i repeat.unmap.fas -s ath -o dir" >ath_NAT_map.sh
qsub -cwd -V -l vf=200M  ath_NAT_map.sh

options:
[mandatory parameters]
	-h|?|--help		help information
	-i reads.fa		Input fasta in miRDeep2 format
	-s species		the abbr of species from "/PUBLIC/source/RNA/smallRNA/version3/6plant_NAT/NAT/PlantNATsDB/spe.list" (a little different from kobas and miRBase)
[optional parameters]
	-o prj			project name, which will create a directory:prj in current dir[default="-s "]
==========================================================================
USAGE
exit 0;
}

my ($help,$reads,$dir,$species);
GetOptions(
	"h|?|help"=>\$help,
	"i=s"=>\$reads,
	"s=s"=>\$species,
	"o=s"=>\$dir,
);


if(defined($help)||!defined($species)||!defined($reads)){
	&usage;
}

my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $bowtie = $Config->{software}->{bowtie1};
my $PlantNATsDB = $Config->{database}->{PlantNATsDB};

my $pwd = `pwd`;
$pwd =~ s/\s+//g;
$reads=get_ful_path($reads);
$dir||=$species;

if(!(-e "$dir")){
	`mkdir -p $dir`;
}
`echo Start Time:`;
`date`;
#chdir("$dir");

open(TAB,"$PlantNATsDB/tab/$species.nat.final.out");

if(!(-e "$dir/input")){
	`mkdir -p $dir/input`;
}
open my $cis,">$dir/input/known_cis_gene.xls";
open my $trans,">$dir/input/known_trans_gene.xls";

print $cis "sense cis gene\tstrand\tantisense cis gene\tstrand\tCIS\tsense start\tsense end\tantisense start\tantisense end\toverlapping length\ttype\n";
print $trans "trans gene one\tstrand\ttrans gene two\tstrand\tTRANS\ttrans gene one start\ttrans gene one end\ttrans gene two start\ttrans gene two end\toverlapping length\n";
my (%cis_gene,%trans_gene);
while(<TAB>){
	chomp;
	my @tmp=split /\t/;
	if($tmp[5] eq "cis"){
		#next if ($tmp[-2]=~/Nearby/i);
		print $cis join("\t",@tmp[1,2,3,4,5,6,7,8,9,10,11]),"\n";
		$cis_gene{$tmp[1]}=1;
		$cis_gene{$tmp[3]}=1;
	}elsif($tmp[5] eq "trans"){
		print $trans join("\t",@tmp[1,2,3,4,5,6,7,8,9,10]),"\n";
		$trans_gene{$tmp[1]}=1;
		$trans_gene{$tmp[3]}=1;
	}
}
close(TAB);
close $cis;close $trans;
`sed 's/\|$species//g' $PlantNATsDB/fa/$species.nat.fas >$dir/input/$species.nat.fas`;
my $in  = new Bio::SeqIO(-file => "$dir/input/$species.nat.fas");
my $out1 = new Bio::SeqIO(-file => ">$dir/input/$species.cis-NAT.fas", -format=>'fasta');
my $out2 = new Bio::SeqIO(-file => ">$dir/input/$species.trans-NAT.fas", -format=>'fasta');
while (my $seq = $in->next_seq()){
	if(defined($cis_gene{$seq->id})){
		$out1->write_seq($seq);
	}
	if(defined($trans_gene{$seq->id})){
		$out2->write_seq($seq);
	}
}

if(!(-e "$dir/output")){
	`mkdir -p $dir/output`;
}
my @types=("cis-NAT","trans-NAT");
my ($rc,$uc);
for my $i(@types){
if(!(-z "$dir/input/$species.$i.fas")){
	`echo bowtie build index for: $species.$i.fas`;
	`$bowtie/bowtie-build $dir/input/$species.$i.fas $dir/input/$species.$i.fas`;
	`$bowtie/bowtie -p 10 -v 0 -k 1 -f $dir/input/$species.$i.fas $reads $dir/output/$i.bwt --un $dir/output/$i.unmap.fas 2>$dir/output/$i.mapping.log`;
	if(!(-z "$dir/output/$i.bwt")){
		`$perlExec $Bin/bwt12collapse.pl $dir/output/$i.bwt >$dir/output/$i.map.fas`;
		`cat $dir/output/$i.map.fas >>$dir/output/NAT.map.fas`;
		`$perlExec $Bin/genebwt12count.pl -i $dir/output/$i.bwt -r $dir/input/$species.$i.fas -t $i -o $dir/output -u -s`;
		`awk \'{if(NR==1){print \"Types\\t\"\$0}else if(NR==2){print \"Mapped reference\\t\"\$0}}\' $dir/output/$i.mapref.stat >$dir/output/$i.map.stat`;
		`awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf(\"Mapped uniq sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"\$i)}printf(\"\\n\");}}\' $dir/output/$i.uc.stat >>$dir/output/$i.map.stat`;
		`awk \'{if(NR==2){for(i=2;i<=NF;i++){total+=\$i}printf(\"Mapped uniq sRNA\\t\"total);for(i=2;i<=NF;i++){printf(\"\\t\"\$i)}printf(\"\\n\");}}\' $dir/output//$i.rc.stat >>$dir/output/$i.map.stat`;
		$reads="$dir/output/$i.unmap.fas";
		if(!defined($rc)){
			($rc,$uc)=("$dir/output/rc.stat","$dir/output/uc.stat");
			`cat $dir/output/$i.uc.stat >$dir/output/uc.stat`;
			`cat $dir/output/$i.rc.stat >$dir/output/rc.stat`;
		}else{
			`awk \'{if(NR>1){print}}\' $dir/output/$i.uc.stat >>$dir/output/uc.stat`;
			`awk \'{if(NR>1){print}}\' $dir/output/$i.rc.stat >>$dir/output/rc.stat`;
		}
	}
}
}
`awk '{if(NR==1){num=NF-1;print}else if(\$1~/+\$/){for(i=1;i<=num;i++){total1[i]+=\$(i+1)}}else if(\$1~/-\$/){for(i=1;i<=num;i++){total2[i]+=\$(i+1)}}else{for(i=1;i<=num;i++){total[i]+=\$(i+1)}}}END{printf("NAT");for(i=1;i<=num;i++){printf("\\t"total[i])}printf("\\n");printf("NAT:+");for(i=1;i<=num;i++){printf("\\t"total1[i])}printf("\\n");printf("NAT:-");for(i=1;i<=num;i++){printf("\\t"total2[i])}printf("\\n");}' $dir/output/uc.stat >$dir/output/NAT.uc.stat`;
`awk '{if(NR==1){num=NF-1;print}else if(\$1~/+\$/){for(i=1;i<=num;i++){total1[i]+=\$(i+1)}}else if(\$1~/-\$/){for(i=1;i<=num;i++){total2[i]+=\$(i+1)}}else{for(i=1;i<=num;i++){total[i]+=\$(i+1)}}}END{printf("NAT");for(i=1;i<=num;i++){printf("\\t"total[i])}printf("\\n");printf("NAT:+");for(i=1;i<=num;i++){printf("\\t"total1[i])}printf("\\n");printf("NAT:-");for(i=1;i<=num;i++){printf("\\t"total2[i])}printf("\\n");}' $dir/output/rc.stat >$dir/output/NAT.rc.stat`;
`cp $reads $dir/output/NAT.unmap.fas`;
`echo End Time:`;
`date`;

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

