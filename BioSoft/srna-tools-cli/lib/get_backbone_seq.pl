#!/usr/bin/perl

###############################################################
#
# This script is part of the srna-tools website.
# It fetches a requested region from one of the genome files
# and writes the sequence into /data/backbone_seq.fasta for
# the SiLoMa tool. 
# 
# Arguments
# --working_dir   Working-directory for this project
# --region        region in the genome in format CHROM:START..END
# --refseq        The path to the selected genome
# --genome_name   Display name of genome, not file name

use strict ;
use warnings;
# On the cluster, this is the path to Bioperl
use lib ("/usr/local/Bioperl/lib/perl5/site_perl/5.8.3", "/usr/lib/perl5/site_perl/5.8.3");
use Bio::DB::Fasta ;
use Bio::SeqIO ;
use Getopt::Long ;

#############################
# Parameters and
# declarations
#############################
my $module_name = "Fetch genomic region" ; # for error log

my $working_dir ;
my $refseq ;
my $region ;
my $update_status_cmd ;
my $data_dir ;
my $output_file = 'backbone_seq.fasta' ;
my $error_file ;
my $genome_name ;

GetOptions(
  'working_dir=s'       => \$working_dir,
	'region=s'            => \$region,
	'refseq=s'            => \$refseq,
  'genome_name=s'       => \$genome_name,
  'update_status_cmd=s' => \$update_status_cmd,
) ;

#############################
# Check parameters and files
#############################

eval {
	die "missing working directory parameter\n" unless $working_dir ;
  die "missing refseq parameter\n" unless $refseq ;
  die "missing region parameter\n" unless $region ;
  
	# Default error file is WORKING_DIR/errors
	$error_file = $working_dir.'/errors' unless $error_file;
  $data_dir = $working_dir.'/data/' ;
  
	# Check files and directories
	die "could not find/read working directory\n" unless -d $working_dir ;
  die "could not find/read genome file\n" unless -r $refseq ;
  die "region paramter not formatted correctly\n" unless $region=~/^.+?:\d+-\d+$/ ;
  die " data directory not found in working directory\n" unless -d $data_dir ;  
} ;
if ($@) {
	handle_fatal("$module_name, error in parameters: $@\n", $error_file) ;
}

# update status on server
eval{ system($update_status_cmd."fetching_region_from_genome") } ;

#############################
# Fetch region from genome
#############################
eval {
  my $db = Bio::DB::Fasta->new($refseq);
  my ($chrom, $start, $end) = ($region=~/^(.+?):(\d+)-(\d+)$/) ;
  my $seq_obj = $db->get_Seq_by_id($chrom);
  die "Could not find chromsome (BAC/scaffold) with ID $chrom in genome $genome_name. If you don't know the exact ID of the required sequence you might want to run the SiLoCo tool first and copy the ID from there." unless $seq_obj ;
  
  my $len = $seq_obj->length ;
  die "The requested region start is outside chromosome $chrom for genome $genome_name (length:$len bp)" if $start > $len ;
  die "The requested region end is outside chromosome $chrom for genome $genome_name (length:$len bp)" if $end > $len ;
  
  my $subseq  = $seq_obj->subseq($start => $end);
  
  die "Could not extract region from chromosome $chrom for genome $genome_name" unless $subseq ;
  
  open (OUT, ">", $data_dir.$output_file) or die "Could not write to output file" ;
  print OUT ">$region ($genome_name)\n$subseq\n" ;
  close OUT ;
} ;
if ($@) {
  handle_fatal("$module_name, error fetching sequence from genome: $@\n", $error_file) ;
}


#############################
# handle_fatal
# This is executed when we catch
# exceptions
# We could add some cleanup code here,
# e.g. close all files and delete samples.
# If the error file can't be opened the script
# will still terminate will exit status > 0
# whih will be picked up by the shell wrapper
# script to gengerate an "errors" file anyway
#############################
sub handle_fatal{
  my $msg = shift ;
  my $error_file = shift ;
  
  open(ERROR_FILE, '>>', $error_file) || die ;
  print ERROR_FILE $msg ;
  close ERROR_FILE ;
  exit 1; # indicate failure
} # handle_fatal
