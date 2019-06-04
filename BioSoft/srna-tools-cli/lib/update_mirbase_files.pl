#!/usr/bin/perl

use strict;
use warnings;
#On the cluster, this is the path to Bioperl
use lib ("/usr/local/Bioperl/lib/perl5/site_perl/5.8.3", "/usr/lib/perl5/site_perl/5.8.3");
use Bio::SeqIO;

########################################################
# This script is part of the srna-tools website.
# It attempts to download the file mature.fa and hairpin.fa
# from mirbase (ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/)
# and saves them as mature.fa and hairpin.fa
# then it generates the files for patman matching with overhang
# by adding "XX" on both sides of each sequence and saves
# these as mature_plusX.fa and hairpin_plusX.fa
########################################################

# warn before continuing
warn <<END ;

###################################
# Warning: this script will download
# new versions of the miRBASE files
# mature.fa and hairpin.fa into the
# current directory and replace existing
# files. It will also generate versions
# called mature_plusX.fa and hairpin_plusX.fa
# for Patman matching with overhangs and
# it generates filtered subsets for all files 
# for Metazoa and Viridiplantae using the
# file organisms.txt from miRBASE.
####################################

CONTINUE [y/n] ?
END

chomp( my $cont = <STDIN>) ;

if ( $cont ne "y" && $cont ne "Y" ) { 
  exit 0 ;
}

# Get the files and uncompress them.
# then make copies for adding the Xs 
# using Perl

my $DOWNLOAD_SUFF=".fa.gz" ;
my $UNCOMPR_SUFF=".fa" ;

my $SANGER_URL="ftp://mirbase.org/pub/mirbase/CURRENT/" ;
my $tax_file = "organism.txt" ;
my $tax_file_rename = "mirbase_organisms.txt" ;
my $cmd ;

# Download the taxonomy files
warn "Fetching $tax_file from $SANGER_URL...\n" ;
$cmd="wget -q ".$SANGER_URL.$tax_file ;
system($cmd) == 0 || die "## Could not download taxonomy file $tax_file - please contact mirbase administrator, it is possible that this file was omitted from the current release\n\n### terminated ###\n" ;
# rename file to mirbase_organisms.txt
warn "done - renaming file to 'mirbase_organisms.txt'\n" ;
rename($tax_file, $tax_file_rename) or die "Could not rename file\n" ;
warn "done\n" ;
warn "\n" ;

# Parse taxonomy file to associate organisms prefixes
# with top level taxonomic unit (Metazoa, Viridiplantae)
my %tax ;
warn "Building table of organisms...\n" ;
open (TAX, '<', $tax_file_rename) or die "Could not open $tax_file_rename\n" ;
while (<TAX>) {
  if (/^(\w+?)\t.+(Metazoa|Viridiplantae)/) {
    $tax{$1} = $2 ;
    warn "$1: $2\n" ;
  }
}
## TODO remvoe as soon as Sam has corrected taxonomy for chlamy  
$tax{'cre'}= 'Viridiplantae' ;
warn "done\n" ;
warn "\n" ;

close TAX ;

# Download the databases
my @dfile_names = ( 'mature','hairpin' );

foreach my $dfile_name (@dfile_names) {
  my $dfile=$dfile_name.$DOWNLOAD_SUFF ;

  warn "Fetching $dfile from $SANGER_URL...\n" ;
  $cmd = "wget -q ".$SANGER_URL.$dfile ;
  system($cmd) == 0 || die "### terminated ###\n" ;
  warn "done\n" ;
  warn "\n" ;
  
  warn "uncompressing file ...\n" ;
  $cmd = "gunzip $dfile" ;
  system($cmd) == 0 || die "### terminated ###\n" ;
  warn "done\n" ;
  warn "\n" ;
  
  my $ufile=$dfile_name.$UNCOMPR_SUFF ; # uncompressed file name
  warn "converting IDs and generating file with XX added to both ends of the sequences ...\n" ;
  
  # shorten fasta IDs to /^>\S/
  # and generate the version with
  # leading/trailing Xs and the sub
  # lists of plant and animal data 
  convert_file($ufile, \%tax) ;
  
  warn "done\n" ;
  warn "\n" ;
}

warn " -- completed successfully --\n" ;
exit 0 ;

########### subroutines

# Generate two files:
# one temp file that contains
# the original sequences but the ID
# lines are shortend to the miRNA IDs
# only (i.e. /^>\S/), this will then
# replace the original input file
# the other one has the sequences with
# "XX" at 5' and 3' end fo Patman matching
# with allowed overlaps of the query
# sequences to the db sequence
sub convert_file {
  my ($infile, $tax_ref) = @_ ; 
  my $tfile="temp_mirbase_conversion_file" ;# a temp file
  
  # file names for:
  # all mirbase (normal and plusX)
  # metazoa
  # viridiplantae
  my ($file_base) = ($infile=~/^(.+?)\./) ;
  my $all_file = $file_base.'_all.fa';
  my $all_fileX =  $file_base.'_all_plusX.fa';
  my $animal_file = $file_base.'_animal.fa';
  my $animal_fileX = $file_base.'_animal_plusX.fa';
  my $plant_file = $file_base.'_plant.fa';
  my $plant_fileX = $file_base.'_plant_plusX.fa';
  
  # open files
  open (ALL_FILE, '>', $all_file) or die "Could not generate file $all_file\n" ;
  open (ALLX_FILE, '>', $all_fileX) or die "Could not generate file $all_fileX\n" ;
  open (ANIMAL_FILE, '>', $animal_file) or die "Could not generate file $animal_file\n" ;
  open (ANIMALX_FILE, '>', $animal_fileX) or die "Could not generate file $animal_fileX\n" ;
  open (PLANT_FILE, '>', $plant_file) or die "Could not generate file $plant_file\n" ;
  open (PLANTX_FILE, '>', $plant_fileX) or die "Could not generate file $plant_fileX\n" ;
   
  my $in = Bio::SeqIO->new( -file => $infile, -format => "Fasta" );
  while (my $is = $in->next_seq) {
    my $seq = $is->seq;
    my $seq_x = "XX".$seq."XX";
    
    my $id = $is->id;
    my $org_prefix ;
    if ($id=~/^(\w+)-/) {
      $org_prefix = $1 ;
    } else {
      die "Could not find organism suffix in ID $id\n#### terminated ####\n"
    }
    
    print ALL_FILE ">$id\n$seq\n";
    print ALLX_FILE ">$id\n$seq_x\n";
    
    if (my $tax = $$tax_ref{$org_prefix}) {
      if( $tax eq 'Metazoa') {
        print ANIMAL_FILE ">$id\n$seq\n";
        print ANIMALX_FILE ">$id\n$seq_x\n";
      } elsif ($tax eq 'Viridiplantae')  {
        print PLANT_FILE ">$id\n$seq\n";
        print PLANTX_FILE ">$id\n$seq_x\n";
      }
    }
    
  }
  unlink $infile or die "Could not remove $infile\n" ;
  close ALL_FILE ;
  close ALLX_FILE ;
  close ANIMAL_FILE ;
  close ANIMALX_FILE ;
  close PLANT_FILE ;
  close PLANTX_FILE ;
}
