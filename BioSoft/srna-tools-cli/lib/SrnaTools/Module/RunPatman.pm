######################################################################
# The UEA sRNA Toolkit: Perl Implementation - A collection of
# open-source stand-alone tools for analysing small RNA sequences.
# Copyright © 2011 University of East Anglia
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
# It runs patman on a given input fasta file in
# WORKING_DIR/data 
# If the number of sequences in the input file
# exceeds a defined maximum (1 million),
# it is split up into chunks with no more 
# than the max number of sequences and
# Patman is run once for each chunk to
# prevent memory problems when running
# Patman with a very large sRNA dataset.
#
# Arguments:
# -w WORKING_DIR
# --patman_bin FULL/PATH/TO/PATMAN/EXECUTABLE
# --refseq FULL/PATH/TO/GENOME/FILE
# -e EDITS                  = -e option for patman (max "edits") defaults to 0
# -s                        = -s option for patman (match plus strand only)
# -o OUTPUTFILE             which will be created in data/
# -i INPUTFILE              assumed to be in data/
# --remove_infiles          removes original sample.infile files
#                           when finished parsing them. This should
#                           only be used when matching to the genome
#                           and not when matching t/rRNAs because the
#                           original files are still required for the
#                           filtering step in the latter case.
# --ensure_result           If set, die if patman result file is empty
#                           This should only be used when matching the genome
#                           and it would be pointless to continue w/o matches 

package SrnaTools::Module::RunPatman ;
use base SrnaTools::Module ;
use strict ;
use warnings;

sub run{
my $self = shift ;

# Modules that may be in non-standard location
require Bio::SeqIO;

#############################
# Parameters and
# declarations
#############################
my $module_name = "Sequence matching" ; # for error log
my $working_dir = $self->job->job_working_dir;
my $data_dir ;
my %srnas ;

# The genome file is usually in the public data dir but
# if we are running in interactive mode on the CLI App then
# it is a complete path. The parameter genome_by_path is set
# in that case
# As a special case, SiLoMa sets the refseq dir explicitly (normally
# to job working_dir/data)
my $refseq_file_name = $self->param('ref_seq');
my $public_data_dir = $self->job->app->path_to('data_dir',$self->job->execution_env);
my $refseq_file ;

if ($self->param('ref_seq_dir')){
  $refseq_file = $self->param('ref_seq_dir').'/'.$refseq_file_name ;
} elsif (substr($refseq_file_name,0,1) eq '/'){
  $refseq_file = $refseq_file_name;
}
else {
  $refseq_file = $refseq_file_name;
  # Try accessing from relative path first, if not, try data dir
  if (-r $refseq_file){
    # do nothing
  } else {
    $refseq_file = $public_data_dir.'/'.$refseq_file_name ;
  }
}

my $outfile_name = $self->param('outfile');
my $infile_name = $self->param('infile');
my $patman_path = 'patman' ; # must be in PATH
my $remove_infiles = $self->param('remove_infiles'); 
my $ensure_result = $self->param('ensure_result'); 
my $use_s = $self->param('same_strand') ;
my $pat_e = $self->param('mismatches') || 0 ; 
my $max_seq_per_file = 1e6 ; # 1million sequences max per patman input file
my $status = $self->param('status_msg') ;

#############################
# Check parameters and files
#############################

eval {
  die "missing working directory parameter\n" unless $working_dir ;
  
  $data_dir = $working_dir.'/data/' ;
  die "missing refseq file parameter\n" unless $refseq_file ;
  die "missing output file name\n" unless $outfile_name ;
  die "missing input file name\n" unless $infile_name ;
  # Check files and directories
  die "could not find/read working directory\n" unless -d $working_dir ;
  die "data directory not found in working directory\n" unless -d $data_dir ;  
  die "could not read fasta file of the sequence you are trying to match sRNAs to. The file may have been removed or renamed by an administrator\n" unless -r $refseq_file ;
  die "patman bin file not accessible or not executable\n" unless $self->binary_in_path('patman');
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error in parameters: $@") ;
}

# update status on server
$self->update_status($status) ;

#############################
# Check size of input file and
# split into smaller chunks if
# it exceeds the maximum number
# of sequences.
#############################
eval{
  my $full_path = $data_dir.'/'.$infile_name ;
  my $pat_s_option = $use_s ? '-s' : '' ; # single strand option for Patman
  my $pat_options = "$pat_s_option -e $pat_e -g 0" ;# patman options
	
  if ($self->file_exceed_max($full_path, $max_seq_per_file)) {
    my $splitup_dir = $self->split_up_infile($data_dir, $infile_name, $max_seq_per_file, $remove_infiles) ;
    
    # Perform PatMan search on each of the smaller files
    opendir(DIR, $splitup_dir) or die "failed to open directory\n";
    while (defined(my $file = readdir(DIR))) {
      next if $file=~/^\./ ;
      my $full_path = $splitup_dir.'/'.$file ;
      my $tmp_outfile = $full_path.'.patout' ;
      $self->do_run_patman($patman_path, $refseq_file, $full_path, $tmp_outfile, $pat_options) ;
    }

    # combine all output files into one in data_dir
    system("cat $splitup_dir/*.patout > $data_dir/$outfile_name") == 0 or die "Could not concatanate result files\n" ;
    system("rm -rf $splitup_dir") ;
    closedir DIR;
  } else { 
    # file does not exceed max num seqs - run patman on whole file
    $self->do_run_patman($patman_path, $refseq_file, $full_path, "$data_dir/$outfile_name", $pat_options) ;
  }

  # At this point we should have a patman result file
  # Terminate if the file doesn't exist.
  # If the option --ensure_result was used, make sure
  # that this file is not zero-length. This is required
  # when we match against the genome and all further
  # steps depend on having matches.
  if (!-e "$data_dir/$outfile_name"){
    die "Patman result file was not created. Please contact administrator.\n"
  }  
  if ($ensure_result && -z "$data_dir/$outfile_name"){
    die "Matching the sRNAs to the genome produced no results. Please check your fasta files (e.g. do they contain adapters, which should have been removed?). Please also make sure that you select the correct genome to match to.\n"
  }
  
} ;         
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while mapping to reference sequence: $@\n") ;
}

return 1;

} #run


######################################################
# Subroutines
######################################################


#############################
# Test if a file contains
# more sequences than the max
# return 1 if it does
#############################
sub file_exceed_max{
  my $self = shift ;
  my $file = shift ;
  my $max_seq_per_file = shift ;
  
  my $cmd = "wc -l $file" ;
  my $output = `$cmd` ;
  my $seq_count ;
  if ($output =~ m/^\s*(\d+)/) {
    $seq_count = $1/2 ;
  } else {
    die "problem retrieving file size. Please contact administrator."
  }
  return ($seq_count > $max_seq_per_file) ? 1 : 0 ;
} # file_exceed_max


#############################
# Split up a large input file 
# into smaller chunks for Patman.
# Create a directory with the
# name ORIGINAL_FILENAME.d
# and write content of orignal
# file to files in that directory
# named ORIGINAL_FILENAME.1 .2 etc
# NOTE: at this stage we assume
# that sequences have been parsed
# before by Bio::SeqIO, so we can
# trust that there are 2 lines in 
# the file for each sequence.
#############################
sub split_up_infile{
  my $self = shift ;
  my $data_dir = shift ;
  my $file = shift ;
  my $max_seq_per_file = shift ;
  my $remove_infiles = shift ;
  
  # create directory for new files
  my $tmp_dir = $data_dir.'/'.$file.'.d' ;
  mkdir $tmp_dir or die "Temp. directory creation failed\n" ;
    
  my $file_num ;
  my $last_file_num = 0 ;
  my $i = 0 ;
  my $out_fh ; # file handle for output files
  
  my $original_file_full_path = $data_dir.'/'.$file ;
  open(FILE,'<',$original_file_full_path) or die "File operation failed\n" ;
  while (my $line=<FILE>) {
    ++$i if substr($line,0,1) eq '>';
    # number of file to write to, starting with 1
    # (will be appended to file name)
    $file_num = int( $i / ($max_seq_per_file + 1) ) + 1;
    if ($file_num != $last_file_num) { # need to start a new file
      if (defined $out_fh) { # close last file
        close $out_fh or die "File operation failed\n" ;
      }
      my $tmp_file = $tmp_dir .'/'.$file.$file_num ;
      open ($out_fh, '>', $tmp_file) or die "Could not create tmp file\n";
      $last_file_num = $file_num ;
    } #if new file_num
    
    print $out_fh $line ;
  } # while <FILE>
  close $out_fh ;
  close FILE ;
  # The original file can be deleted if
  # "remove_infiles" option selected
  unlink $original_file_full_path if $remove_infiles ;
  
  # return the name of the dir
  # containing the splitup fasta files
  return $tmp_dir ;

} #split_up_infile


#############################
# Run patman with given options
#############################
sub do_run_patman{
  my $self = shift ;
  my ($patman_path, $refseq_file, $full_path, $outfile, $pat_options) = @_ ;
  my $cmd ="$patman_path -D $refseq_file -P $full_path $pat_options -o $outfile" ;
  my $output = `$cmd` ;
  die "Problem running PatMaN\n" if $output ;
}

1;