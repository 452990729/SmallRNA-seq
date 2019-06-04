#!/usr/bin/perl -w

use Getopt::Long;
use strict;
use Config::Tiny;


my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../BioModule/Settings/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $pythonExec = $Config->{srnaenv}->{python_v276};
my $kobas_v20120208 = $Config->{database}->{kobas_v20120208};
my $blast = $Config->{software}->{blast2};
my $blast1 = $Config->{software}->{blast1};
my $help;
my $input;
#my $blast_sh;
#my $annot_sh;
#my $identify_sh;
my $output;
my $sample;
my $species = "ko";
my $database = "r/B/g/f/K/o/N/p/R/n/b/k";
my $blast_program = "$blast/blastx";
my $e_value = 1e-5;
my $method_for_identify = "h";
my $false_discovery_rate = "QVALUE";
my $number_of_processors = 1;
my $standard = "not standard";
my $sp_list_dir = "$kobas_v20120208/sp.list";
my $sp_data_dir = "$kobas_v20120208/seq_pep_v2/";
my $dirBlastResult;
my $dirAnnotateResult;
my $dirIndentifyResult;
my $dirAnnotate = "$Bin/annotate.py";
my $dirIdentify = "$Bin/identify.py";

GetOptions("h|help" => \$help,
#	    "bsh=s" => \$blast_sh,
#	    "ash=s" => \$annot_sh,
#	    "ish=s" => \$identify_sh,
           "i=s" => \$input,
           "o=s" => \$output,
           "s=s" => \$species,
           "sample-name=s" => \$sample,
           "e=f" => \$e_value,
           "m=s" => \$method_for_identify,
           "f=s" => \$false_discovery_rate,
           "n=i" => \$number_of_processors,
           "d=s" => \$database,
           "blast-result=s" => \$dirBlastResult,
           "annotate-result=s" => \$dirAnnotateResult,
           "bp=s" => \$blast_program,
           "dir-species-list=s" => \$sp_list_dir,
           "dir-species-data=s" => \$sp_data_dir) or die "Cannot get the inputs: $!";

my $sp_name_list;
my $sp_name_data;
#my $sp_name_data_dir; by lixr

if (!defined($output) || !defined($input) || defined($help))
{
	&usage;
	exit 0;
}

my $dirTemp;
$dirTemp = "$output/temp";
$dirTemp =~ s/\/\//\//;
unless (-d $dirTemp)
{
        !system "mkdir -p $dirTemp" or warn "Something went wrong!";
}

open SP_LIST, $sp_list_dir or die "Cannot open the file: $!";
while(<SP_LIST>)
{
  chomp;
  if(/\s+($species)$/)
  {
    #s/\s*($species)$//;
    s/\s+($species)$//;
    $sp_name_list = $_;
    last;
  }
}
close SP_LIST;

opendir SP_DATA, $sp_data_dir or die "Cannot open the folder: $!";
foreach my $file (readdir SP_DATA)
{
  if ($file =~ /^($sp_name_list).+?fasta$/)
  {
    $sp_name_data = $file;
    last;
  }
}
closedir SP_DATA;

$sp_data_dir = "$sp_data_dir\/$sp_name_data";
print "$sp_data_dir\n";

$dirBlastResult = "$dirTemp/$sample.blast";
$dirAnnotateResult = "$dirTemp/$sample.annotate";
$dirIndentifyResult = "$output/$sample.indentify.xls";
$dirIndentifyResult =~ s/\/\//\//;

open TB, ">$output/blast.sh";
open TA, ">$output/annot.sh";
open TI, ">$output/identify.sh";
print TB "$blast1/blastall -p $blast_program -d $sp_data_dir -a $number_of_processors -e $e_value -i $input -o $dirBlastResult -m 7\n";

if ($species eq "000" || $species eq "001" || $species eq "002")
{
	print TA "$pythonExec $dirAnnotate -i $dirBlastResult -t blastout:xml -s ko -n $number_of_processors -o $dirAnnotateResult\n";
}
else
{
	print TA "$pythonExec $dirAnnotate -i $dirBlastResult -t blastout:xml -s $species -n $number_of_processors -o $dirAnnotateResult\n";
}

if ($species eq "000" || $species eq "001" || $species eq "002")
{
	print TI "$pythonExec $dirIdentify -f $dirAnnotateResult -d K -m $method_for_identify -n $false_discovery_rate -b ko -s $standard -o $dirIndentifyResult\n";
}
else
{
	print TI "$pythonExec $dirIdentify -f $dirAnnotateResult -d K -m $method_for_identify -n $false_discovery_rate -b $species -s $standard -o $dirIndentifyResult\n";
}

#------------------------------------------------

sub usage
{
	print STDERR <<USAGE;
===============================================================================
kobas2.0
Function: KOBAS 2.0 is an update of KOBAS (KEGG Orthology Based Annotation System). Its purpose is to identify statistically enriched related pathways and diseases for a set of genes or proteins, using pathway and disease knowledge from
multiple commonly used databases.
-------------------------------------------------
Usage: runKOBAS.pl [options]
options:
	-i	<string>	input .fa files
	-o	<string>	directory of output files
	-s	<string>	abbr for species
	--sample-name	<string>	the sample name you analyze
	-e	<float>		e value
	-m	<string>	statistic method for identify, b is binomial
				test, c is chisquare test, f is fisher exact
				test, h is hypergeometric test and x is
				frequency list. Default: hypergeometric test.
	-f	<string>	false discovery rate correction method: QVALUE,
				BH, BY or none. Default: QVALUE.
	-n	<int>		number of processors
	-d	<string>	selected databases, one letter abbreviation
				separated by "/", K is KEGG PATHWAY, n is PID
				Curated, b is PID BioCarta, r is PID Reactome,
				B is BioCyc, R is Reactome, p is PANTHER, k is
				KEGG DISEASE, f is FunDO, o is OMIM, g is GAD, N
				is NHGRI, default r/B/g/f/K/o/N/p/R/n/b/k
	--dir-species-list	<string>	species list
	--dir-species-data	<string>	directory of KEGG database
	-p	<string>	the blast program you choose for blastall.

				  blastn: compares a nucleotide query sequence
				against a nucleotide sequence database;
				  blastp: compares a protein query sequence
				against a protein sequence database;
				  blastx: compares a nucleotide query sequence
				translated in all reading frames against a
				protein sequence database. You could use this
				option to find potential translation products
				of an unknown nucleotide sequence;
				  tblastn: compares a protein query sequence
				against a nucleotide sequence database
				dynamically translated in all reading frames;
				  tblastx: compares the six-frame translations
				of a nucleotide query sequence against the
				six-frame translations of a nucleotide sequence
				database;
				Note: gapped alignments are not available with
				tblastx
				  psiblast: iterative form of blastp in which a
				profile is created from the amino acid query and
				nth set of results (meeting the Psi-Expectation)
				and resubmitted
USAGE
}
