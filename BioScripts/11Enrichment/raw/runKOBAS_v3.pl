#!/usr/bin/perl -w

use Getopt::Long;
use strict;
use FindBin '$Bin';
use Config::Tiny;


my $Config = Config::Tiny -> new();
$Config = Config::Tiny ->read("$Bin/../../Pipeline/config.ini");
my $perlExec = $Config->{srnaenv}->{perl_v5182};
my $pythonExec = $Config->{srnaenv}->{python_v276};
my $R_v2153 = $Config->{srnaenv}->{R_v2153};
my $kobas_v20140801 = $Config->{database}->{kobas_v20140801};
my $LIBRARY_PATH = $Config->{srnaenv}->{LIBRARY_PATH};


my $help;
my $input;
my $output;
my $sample;
my $species = "ko";
my $e_value = 1e-5;
my $method_for_identify = "h";
my $false_discovery_rate = "QVALUE";
my $number_of_processors = 1;
my $standard = "not standard";
my $sp_data_dir = "$kobas_v20140801/seq_pep/";
my $BlastResult;
my $dirAnnotateResult;
my $dirIndentifyResult;
my $dirAnnotate = "$kobas_v20140801/scripts/annotate.py";
my $dirIdentify = "$kobas_v20140801/scripts/identify.py";
my $split = "$Bin/split_blastout_v2.pl";
GetOptions("h|help" => \$help,
           "i=s" => \$input,
           "o=s" => \$output,
           "s=s" => \$species,
           "sample-name=s" => \$sample,
           "e=f" => \$e_value,
           "m=s" => \$method_for_identify,
           "f=s" => \$false_discovery_rate,
           "n=i" => \$number_of_processors,
           "blast-result=s" => \$BlastResult,
#           "dir-species-list=s" => \$sp_list_dir,
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

#open SP_LIST, $sp_list_dir or die "Cannot open the file: $!";
#while(<SP_LIST>)
#{
#  chomp;
#  if(/\s+($species)$/)
#  {
#    #s/\s*($species)$//;
#    s/\s+($species)$//;
#    $sp_name_list = $_;
#    last;
#  }
#}
#close SP_LIST;

opendir SP_DATA, $sp_data_dir or die "Cannot open the folder: $!";
foreach my $file (readdir SP_DATA)
{
  if ($file =~ /^($species).+?fasta$/)
  {
    $sp_name_data = $file;
    last;
  }
}
closedir SP_DATA;

$sp_data_dir = "$sp_data_dir\/$sp_name_data";
print "$sp_data_dir\n";

my $dirBlastResult = "$dirTemp/$sample.blast";
$dirAnnotateResult = "$dirTemp/$sample.annotate";
$dirIndentifyResult = "$output/$sample.identify.xls";
$dirIndentifyResult =~ s/\/\//\//;

open TA, ">$output/annot.sh";
open TI, ">$output/identify.sh";
print TA "$perlExec $split -xml $BlastResult -geneid $input -out $dirBlastResult\n";
print TA "export PATH=$pythonExec:$R_v2153:\$PATH\n";
print TA "export PYTHONPATH=$kobas_v20140801/src:\$PYTHONPATH\n";
print TA "export LD_LIBRARY_PATH=$LIBRARY_PATH/R_v2153:\$LD_LIBRARY_PATH\n";
if ($species eq "000" || $species eq "001" || $species eq "002")
{
	print TA "$pythonExec $dirAnnotate -i $dirBlastResult -t blastout:xml -s ko -n $number_of_processors -o $dirAnnotateResult\n";
}
else
{
	print TA "$pythonExec $dirAnnotate -i $dirBlastResult -t blastout:xml -s $species -n $number_of_processors -o $dirAnnotateResult\n";
}

print TI "export PATH=$pythonExec:$R_v2153:\$PATH\n";
print TI "export PYTHONPATH=$kobas_v20140801/src:\$PYTHONPATH\n";
print TI "export LD_LIBRARY_PATH=$LIBRARY_PATH/R_v2153:\$LD_LIBRARY_PATH\n";
if ($species eq "000" || $species eq "001" || $species eq "002")
{
	print TI "$pythonExec $dirIdentify -f $dirAnnotateResult -d K -m $method_for_identify -n $false_discovery_rate -b ko -s $standard -o $dirIndentifyResult\n";
}
else
{
	print TI "$pythonExec $dirIdentify -f $dirAnnotateResult -d K -m $method_for_identify -n $false_discovery_rate -b $species -o $dirIndentifyResult\n";
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
	-i	<string>	input diffgene.id files
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
	--dir-species-list	<string>	species list
	--dir-species-data	<string>	directory of KEGG database
-------------------------------------------------
USAGE
}
