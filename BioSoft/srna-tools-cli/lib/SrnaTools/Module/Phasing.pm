######################################################################
# The UEA sRNA Toolkit: Perl Implementation - A collection of
# open-source stand-alone tools for analysing small RNA sequences.
# Copyright � 2011 University of East Anglia
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
######################################################################


#####################################################################
#
# This is a module of the sRNA tool kit.
# It runs the ta-siRNA analysis
# TODO This class has been converted from a standalone script and might
# still need some cleanup (e.g. passing of parameters)
#####################################################################

package SrnaTools::Module::Phasing ;
use base SrnaTools::Module ;
use strict;
use warnings;

sub run{
my $self = shift ;

my $input = $self->param('input_file');
my $cutoff = $self->param('pval');
my $got_results = 0;
my $min_abundance = $self->param('min_abundance');
my $working_dir = $self->job->job_working_dir ;

open (LOCUSLIST, ">$working_dir/results/locuslist.csv") or SrnaTools::Exception::ModuleExecution->throw("Couldn't create output file for predicted TAS locus list\n") ;
open (SRNAS, ">$working_dir/results/srnas.txt") or SrnaTools::Exception::ModuleExecution->throw("Couldn't create output file for predicted ta-siRNA list\n") ;

if (!$cutoff){
  $cutoff = 0.01;
}
my $found = 0;
unless ($input){
  print STDERR "No input file given - exiting\n";
  exit();
}


my $last_start = 0;
my $last_chr = "";
my %seen;
my $output;
my $maxp=100;
my %sh;
my %seq_info;
my %seq_abundance;
my %seq_output;
open(DATA, "$input") or SrnaTools::Exception::ModuleExecution->throw("Couldn't open patman output file\n") ;
print LOCUSLIST "Using cutoff p-value of $cutoff\n\n";
print LOCUSLIST "Chromosome,Start postition,End position,# sequences,# phased sequences,p-value\n";
print SRNAS "List of phased sRNAs in each of the predicted loci\n";
while(<DATA>){
	if ($_=~m/^(\S+)\s*.*\t(\w+)\_\S+\((\d+)\)\t(\d+)\t(\d+)\t(\S+)\t\d+/){
	my $chr = $1;
	my $cor = $4;
	my $str = $6;
	my $seq = $2;
	my $end = $5;
	my $abundance = $3;
	chomp ($chr,$cor,$str,$seq,$abundance);

		if ($str eq '+'){
		$str = "1";
		}
		else{
		$str = "-1";
		$cor = $end;
		}
		if (length($seq) != 21){
		next;
		}
		else{
		my $line = "$chr\t$cor\t$str\t$seq";
		my $sp=$chr.",".$cor.",".$str;
			if ($abundance >= $min_abundance){
			$sh{$sp}=$line;
			$seq_info{$sp}=$seq;
			$seq_abundance{$sp}=$abundance;
			}
		}
	}
}
close DATA;
my $seq_ref = 0;
foreach my $ccs (sort keys %sh){
chomp $ccs;
my $n=0;
my $k=0;
my ($chr, $cor, $str)=split(/,/, $ccs);
my ($abundance,$sequence, $outline);
if ($str==1){ # small RNAs on the Watson strand
my $ss=$cor;
my $ee=$cor+230;
	for (my $i=$ss; $i<=$ee; $i++){ # calculate n on the sense strand
	my $x=$chr.",".$i.","."1";
		if (exists $sh{$x}){
		$n=$n+1;
		}
	}
	for (my $i=$ss; $i<=$ee; $i=$i+21){ # calculate k on the sense strand
	my $x=$chr.",".$i.","."1";
		if (exists $sh{$x}){
		$k=$k+1;
		$abundance = $seq_abundance{$x};
		$sequence = $seq_info{$x};
		$outline = "$sequence($abundance)\t$x";
		push (@{$seq_output{$ccs}}, $outline);
		}

	}
	for (my $j=$ss-2; $j<=$ee-2; $j++){ # calculate n on the antisense strand
	my $x=$chr.",".$j.","."-1";
		if (exists $sh{$x}){
		$n=$n+1;
		}

	}
	for (my $j=$ss+18; $j<=$ee-2; $j=$j+21){ # calculate k on the antisense strand
	my $x=$chr.",".$j.","."-1";
		if (exists $sh{$x}){
		$k=$k+1;
		$abundance = $seq_abundance{$x};
		$sequence = $seq_info{$x};
		$outline = "$sequence($abundance)\t$x";
		push (@{$seq_output{$ccs}}, $outline);
		}
	}

}
elsif ($str==-1){ # small RNAs on the Crick strand
my $ee=$cor;
my $ss=$cor-230;
	for (my $i=$ss; $i<=$ee; $i++){ # calculate n on the sense strand
	my $x=$chr.",".$i.","."-1";
		if (exists $sh{$x}){
		$n=$n+1;
		}
	}
	for (my $i=$ss+20; $i<=$ee; $i=$i+21){ # calculate k on the sense strand
	my $x=$chr.",".$i.","."-1";
		if (exists $sh{$x}){
		$k=$k+1;
		$abundance = $seq_abundance{$x};
		$sequence = $seq_info{$x};
		$outline = "$sequence($abundance)\t$x";
		push (@{$seq_output{$ccs}}, $outline);
		}	
	}
	for (my $j=$ss+2; $j<=$ee+2; $j++){ # calculate n on the antisense strand
	my $x=$chr.",".$j.","."1";
		if (exists $sh{$x}){
		$n=$n+1;
		}

	}
	for (my $j=$ss+2; $j<=$ee+2; $j=$j+21){ # calculate k on the antisense strand
	my $x=$chr.",".$j.","."1";
		if (exists $sh{$x}){
		$k=$k+1;
		$abundance = $seq_abundance{$x};
		$sequence = $seq_info{$x};
		$outline = "$sequence($abundance)\t$x";
		push (@{$seq_output{$ccs}}, $outline);
		}
	}

}

my $p=0;

for (my $w=$k; $w<=21; $w++){ # calculate p-value from n and k
my $c=1;
my $rr=1;
my $rw=1;
	for (my $j=0; $j<=$w-1; $j++){
	$c=$c*($n-$j)/($j+1);
	}
	for (my $x=0; $x<=$w-1; $x++){
	$rr=$rr*(21-$x)/(461-$x);
	}
	for (my $y=0; $y<=$n-$w-1; $y++){
	$rw=$rw*(440-$y)/(461-$w-$y);
	}
my $pr=$c*$rr*$rw;
$p=$p+$pr;
}

$p = sprintf("%e", $p);

if ($p <= $cutoff){
my ($chr, $start, $str, $seq)=split(/\t/, $sh{$ccs});
	if ( (abs($start-$last_start) > 251) || ($chr ne $last_chr) ){
		if ($output){
		print LOCUSLIST "$output";
		print SRNAS "$output\n";
			foreach my $el (@{$seq_output{$seq_ref}}){
			print SRNAS "$el\n";
			}
		print SRNAS "\n";
		}
		if ($n > 2){
		my $end;
		my $st = $start;
			if ($str == -1){
			$end = $start+21;
			$st = $start-230;
			}
			else{
			$end = 251+$start;
			}
		$output = "$chr,$st,$end,$n,$k,$p\n";
		$seq_ref = $ccs;
		}
		else{
		undef $output;
		}
	$maxp = $p;
	}
	else{
		if (($p < $maxp) && ($n > 2)){
		my $end;
		my $st = $start;
			if ($str == -1){
			$end = $start+21;
			$st = $start-230;
			}
			else{
			$end = 251+$start;
			}
		$output = "$chr,$st,$end,$n,$k,$p\n";
		$seq_ref = $ccs;
		$maxp = $p;
		}
	}
$last_start = $start;
$last_chr = $chr;
$got_results++;
}
}
if ($output){
print LOCUSLIST "$output";
print SRNAS "$output\n";
foreach my $el (@{$seq_output{$seq_ref}}){
  print SRNAS "$el\n";
}
print SRNAS "\n";
$seq_ref = 0;
}

if (!$got_results){
print LOCUSLIST"No results found at using this cutoff (please make sure you are using a redundant sequence set - all reads of abundance 1 are ignored)\n";
print SRNAS "No results found at using this cutoff (please make sure you are using a redundant sequence set - all reads of abundance 1 are ignored)\n";
}
close SRNAS;
close LOCUSLIST;

return 1;
} # run

1;