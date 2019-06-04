use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#BY JC
my ($RNAhybrid,$PITA,$miRanda);

GetOptions(
	"RNAhybrid=s" =>\$RNAhybrid,
	"PITA=s" =>\$PITA,
	"miRanda=s" =>\$miRanda,
	);

my $usage="
USAGE:
	perl $0 
		-RNAhybrid	RNAhybrid result
		-PITA		PITA result
		-miRanda	miRanda result\n";
for ( ($RNAhybrid,$PITA,$miRanda) ){
	die $usage unless defined $_;
}
#get RNAhybird result
my %all_target;
my %RNAhybrid;
sub RNAhybrid{
	my ($target,$miRNA);
	open IN ,$RNAhybrid;
	while(<IN>){
		chomp;
		if (/^target:\s*(.*)$/){
			$target=$1;
			while (<IN>){
				chomp;
				if (/^miRNA :\s*(.*)$/){
					$miRNA=$1;
					$RNAhybrid{"$target::$miRNA"}=1;
					$all_target{"$target::$miRNA"}=1;
					last;
					}
				}
		}
	}
	close IN;
#	print Dumper \%RNAhybrid;
}

#get PITA result
my %PITA;
sub PITA{
	open IN,$PITA;
	<IN>;
	while (<IN>){
		chomp;
		my $miRNA=(split /\t/)[1];
		my $target=(split /\t/)[0];
		$PITA{"$target::$miRNA"}=1;
		$all_target{"$target::$miRNA"}=1;
	}
	close IN;
#	print Dumper \%PITA;
}

#get miRanda result
my %miRanda;
sub miRanda{
	open IN,$miRanda;
	while (<IN>){
		chomp;
		my @temp=split /\t/;
		$miRanda{"$temp[1]::$temp[0]"}=1;
		$all_target{"$temp[1]::$temp[0]"}=1;
	}
	close IN;
}

sub common_target {
	open OUT, ">commom_target.xls";
	open OUT_2,">all_target.xls";
	print OUT "Target::miRNA\tRNAhybrid\tPITA\tmiRanda\n";
	RNAhybrid();
	print "RNAhybrid finished\n";
	miRanda();
	print "miRanda finished\n";
	PITA();
	print "PITA finished\n";
	for (sort keys %all_target){
		my $result="$_\t";
		$result.= exists $RNAhybrid{$_} ? "1\t" : "0\t";
		$result.= exists $PITA{$_} ? "1\t" : "0\t";
		$result.= exists $miRanda{$_} ? "1\n" : "0\n";
		print OUT "$result";
		if ( exists $RNAhybrid{$_} and exists $PITA{$_} and exists $miRanda{$_}  ){
			my ($target,$miRNA)=$_=~ /(\S+)::(\S+)/;
			die "Can not find target and miRNA pairs " unless $target and $miRNA;
			print OUT_2 "$miRNA\t$target\n";
		}
	}

}
common_target();
