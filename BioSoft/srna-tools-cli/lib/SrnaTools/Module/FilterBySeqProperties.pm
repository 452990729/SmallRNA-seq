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
# It extracts sequences from input files and (if any)
# public datasets requested. Filter sequences
# to exclude those that are
# - too short/long
# - low complexity (less than 3 different bases)
# - contain other than AGCT (translate U to T)
# not valid fasta format
#
# NOTE This has been converted from a non-OO version
# of the script and still needs some cleaning up
#
# Old Documentation
# All sequences from all files will be written
# to one new file in non-redundant format
# and counts for each sample will be stored in the
# ID lines in the format:
#>SEQUENCE_SAMPLENAME(ABUNDANCE)[_SAMPLENAME(ABUNDANCE)]
#
# -o OUTPUTFILE  Name of the output file (non-redundant combined fasta list)
# --lc           Apply low-complexity filter
# --min_length   minimum sRNA length
# --max_length   maximum sRNA length
# --max_num_seq  terminate program if the number of valid unique sequences
#                exceeds this number. If this option is not used, the program
#                will not terminate itself.
# -p FILE        Name(s) of GEO files in publid dir to be included in analysis
# --public_dir   Path to public (GEO) sample files
# --filter_name  how to call this filter in config file (will appear later in 
#                stats in tools that make use of this feature. Defaults to 'filter valid sequence'

package SrnaTools::Module::FilterBySeqProperties ;
use base SrnaTools::Module ;
use strict;
use warnings;

sub run{
my $self = shift ;

# Modules that may be in non-standard location
require Bio::SeqIO;

#############################
# Declarations
#############################
my $module_name = "Filtering by sequence properties" ; # for error log
my $filter_low_comp  = $self->param('filter_low_comp');

# The sRNA size range
# use some generous defaults in case they
# were not set
my $min_length = $self->param('min_length') || 1;
my $max_length = $self->param('max_length') || 40;

# Use this to let files pass where no sequence 
# occurrs more than once (off by default)
my $allow_nr = $self->param('allow_nr') ;

my $max_num_seq = $self->job->max_number_uniq_seqs; 

my $working_dir = $self->job->job_working_dir;  
my $outfile = $self->param('outfile') ;
my $data_dir ;
my $config_file  ;
my %in_stats ;
my %valid_stats ;
# how to call this filtering step in config file
my $config_filter_name='filter valid sequence' ; 

my %filter_args ;

#############################
# Check parameters and files
#############################

eval {
  die " missing working directory parameter\n" unless $working_dir ;
  die "missing output file name\n" unless $outfile ;
  
  $data_dir = $working_dir.'/data/' ;
  
#   if (!$public_dir && @public_file_names) {
#     die " missing public data directory parameter while using public data\n" ;
#   }
  die " could not find/read working directory\n" unless -d $working_dir ;
  die " data directory not found in working directory\n" unless -d $data_dir ;  
#   die " could not find/read public data directory\n" if $public_dir && ! -d $public_dir ;
  
  # TODO, get rid of this
  %filter_args = ( lc          => $filter_low_comp, 
                   min_length  => $min_length, 
                   max_length  => $max_length,
                   max_num_seq => $max_num_seq );
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error in parameters: $@") ;
}

# update status on server
$self->update_status("extracting_valid_sequences") ;

#############################
# Filter sequences 
# from user files and public
# datasets and convert into
# a single nr list with
# counts for each sample.
#############################

my %srnas ; # nr list of sequences with counts
my %invalid_srnas; # nr list of invalid sequences

# samples will have a short name "S1,S2..." in the nr file
# to save space, here we keep the association between this short
# name and the name of the file parameter (e.g. srna-file)
# and put in the cache for other modules. NOTE This is NOT the
# name of the uploaded file on the user's disk!
my %sample_to_filenames ;

eval{
  # Filter user files
  # Traverse list of data file parameters from tool
  # config and process the file, if present, that
  # corresponds to that. e.g.: many tools will have a 
  # single 'srna_file', so we look for srna_file.uncompressed
  my $sample_count = 1 ;
  foreach my $file_basename (@{$self->job->data_file_params}) {
    my $infile = $data_dir.'/'.$file_basename.'.uncompressed' ;
    my $sample_name = 'S'.$sample_count ;
    $sample_to_filenames{$sample_name} = $file_basename;
    ++$sample_count;
    next unless -e $infile ;
    
    # Filter the sequences and add to non-redundant
    # list with counts for each sample
    # Return total counts (nr count is 
    # from valid sequences only)
    # TODO in urgent need of cleaning up (make use of instance vars)
    my ($total_seq, $nr_seq, $total_seq_valid, $nr_seq_valid, $invalid_char, $invalid_low_comp, $invalid_size) = $self->filter_sequences($infile, $sample_name, \%srnas, \%invalid_srnas, \%filter_args);
    unless ($total_seq) { # we have no sequences at all
      if (-z $infile) {
        die "No sequences could be extracted from $file_basename - file is empty. Please try to submit the job again and check your sequence files carefully if you are getting the same error message. Please contact us if you are sure that your file is ok.\n" ; 
      } else {
        die "No sequences could be extracted from file $file_basename (file size: ".((-s $infile)/1000)." kb). Please check that all uploaded sRNA sequence files are in valid FASTA format.\n" ; 
      }
    }
    # we do have sequences but do we have any valid ones? 
    unless ($total_seq_valid) {
      die "No valid sequences could be extracted from file $file_basename. Number of invalid sequences: $invalid_char contained characters other than A,G,C or T (sequences expected in DNA form), $invalid_low_comp contained less than 3 different bases (low-complexity sequences) and $invalid_size were < ".$filter_args{min_length}."nt or > ".$filter_args{max_length}."nt long. Number of sequences in total: $total_seq. Please make sure you upload PLAIN TEXT FASTA files only (no Word or rich-text files). It is best to use a plain text editor, such as Notepad on Windows or TextEdit on MAC OS X to edit files before uploading.\n" ; 
    }
    
    # Uploaded sequences should always be
    # in redundant format so that we can
    # obtain counts. Terminate if nr count
    # and red. count are equal, which 
    # should be impossible for a redundant
    # 454 or SOLEXA dataset
    if ($total_seq_valid == $nr_seq_valid and !$allow_nr ) {
      die "The total number of valid reads and the number of non-redundant valid reads are equal in file $file_basename. This is probably due to your input set having been made non-redundant before uploading. The tools on this site need to obtain counts for each sRNA sequence, so it is important that you don't make your sequences non-redundant before uploading.\n"  ;
    }
    
    # add to stats
    $in_stats{$sample_name} = {total => $total_seq, nr => $nr_seq } ;
    $valid_stats{$sample_name} = {total => $total_seq_valid, nr => $nr_seq_valid } ;
    
    # We don't need the original file anymore
    system("rm -f $infile") ;
  }
	
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while filtering sequences: $@") ;
}


#############################
# Write stats to config file
#############################

eval{
  $self->write_stats_to_config(\%in_stats, \%valid_stats, $config_filter_name, \%sample_to_filenames) ;
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while writing to config file: $@") ;
}


#############################
# Write nr list to new file
# nr_seqs.infile in data dir
#############################
eval{
  $self->write_nr_seq_file(\%srnas, $data_dir, $outfile) ;
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while writing non-redundant sequence list to tmp file: $@") ;
}

return 1 ;

} # run



############################# 
# Filter sequences
# filter out if outside
# accepted range and (optional)
# low complexity.
# Store valid seqs in nr list
# with abundance counts for
# current sample.
#############################
sub filter_sequences{
  my $self = shift ;
  my $infile = shift;
  my $sample_name = shift ;
  my $srnas_ref = shift ;
  my $invalid_srnas_ref = shift;
  my $filter_args_ref = shift ;
  
  # keep counts
  my $valid_seq = 0 ;
  my $total_seq = 0 ;
  my $total_nr_seq = 0 ;
  my $valid_nr_seq = 0 ;
  my $invalid_char = 0 ;
  my $invalid_low_comp = 0 ;
  my $invalid_size = 0 ;
  my $invalid_nr_seq = 0 ;

  die "Sequence input file $infile not readable - please contact administrator\n" unless -r $infile ;
  my $in_seqio = Bio::SeqIO->new ( -format => 'fasta', -file => $infile ) || die " problem creating Bio::SeqIO object\n";
    
  while (my $inseq = $in_seqio->next_seq) {
    my $seq = uc $inseq->seq ;
    $seq=~tr/U/T/ ;
    my $size = $inseq->length;
    ++$total_seq ;
    
    # If this sequence exists already
    # in the hash of valid nr seqs
    # just update the counters.
    # If this sequence has no count
    # for the current sample, increase
    # the nr-seq count for this sample
    if (defined $$srnas_ref{$seq} ) {
      ++$valid_nr_seq unless defined $$srnas_ref{$seq}{$sample_name} ;
      ++$valid_seq ; 
      ++$$srnas_ref{$seq}{$sample_name} ;
      next ;
    } elsif ( defined $$invalid_srnas_ref{$seq} ) {
      next ; # we have already rejected this one
    }
    
    # Sequence is not in the list yet:
    # check whether it is valid and
    # add to list (with sample abundance
    # count) if it is.
    my $is_invalid = 0;
    if ($seq!~/^[AGCT]+$/) {  # only allow AGC and T
      ++$invalid_char ;
      $is_invalid = 1;
    } elsif ($self->is_low_complexity($seq) && $$filter_args_ref{lc} ) { # test for only one or two different bases
      ++$invalid_low_comp ;
      $is_invalid = 1 ;
    } elsif ( $size < $$filter_args_ref{min_length} ||  # reject incorrect size
              $size > $$filter_args_ref{max_length} ) {
      ++$invalid_size ;
      $is_invalid = 1 ;
    }
    
    if ($is_invalid) {
      $$invalid_srnas_ref{$seq} = 1;
      ++$invalid_nr_seq ;
      next ;
    }
    
    # If we are here, sequence is valid
    # and not in the list yet, so it
    # is a new valid "nr" sequence
    ++$valid_nr_seq ;
    ++$valid_seq;
    ++$$srnas_ref{$seq}{$sample_name} ;
    
    # To avoid memory problems, terinate the program
    # if a given number of sequences is reached
    if (defined $$filter_args_ref{max_num_seq} && $valid_nr_seq > $$filter_args_ref{max_num_seq}) {
      die "The maximum number of unique sequences ($$filter_args_ref{max_num_seq} was exceeded by the selected/uploaded samples. Due to memory limitations, we are unable to process datasets that exceed this maximum. We apologise for the inconvenience." ;
    }
  } # parser
  
  # For now return empty 
  $total_nr_seq = $invalid_nr_seq + $valid_nr_seq ; 
  
  return ($total_seq, $total_nr_seq, $valid_seq, $valid_nr_seq, $invalid_char, $invalid_low_comp, $invalid_size) ;

} # upload_file


#############################
# Test a sequence for
# being low complexity
# i.e. less than 3 different
# bases to exclude things like
# AGAGAGAGAGAG or AAAAAAAAAAA
# return 1 if it is low complexity
# zero if it isn't
sub is_low_complexity{
  my $self = shift ;
  my $seq = shift ;
  
  my %bases ;
  foreach (split('',$seq)) {
    $bases{$_} = 1 ;
  } 
  
  my $bcount = keys %bases;
  if ($bcount <3) { # low complexity sequence
    return 1 ;
  } else { # not low complexity
    return 0;
  }
} # is_low_complexity


#############################
# This used to be in a config
# file but now it is held in 
# the cache managed by the job
#############################
sub write_stats_to_config{
  my $self = shift ;
  my $in_stats_ref = shift ;
  my $valid_stats_ref = shift ;
  my $filter_name = shift ;
  my $sample_to_filenames_ref = shift ;
 
  push @{$self->job->cache->{read_counts}}, {input => $in_stats_ref} ;
  push @{$self->job->cache->{read_counts}}, {$filter_name => $valid_stats_ref} ;
  $self->job->cache->{sample_to_filenames} = $sample_to_filenames_ref ;
 
} # write_stats_to_config


#############################
# write new nr_seq file
# nr_seqs.infile in data dir
#############################
sub write_nr_seq_file{
  my $self = shift ;
  my $srnas_ref = shift ;
  my $data_dir = shift ;
  my $outfile = shift ;
  
  my $file = $data_dir.'/'.$outfile ;
  open(FILE, '>', $file) or die "Could not create file" ;
  
  foreach my $seq (keys %$srnas_ref) {
    print FILE ">$seq" ;
    foreach my $sample (keys %{$$srnas_ref{$seq}} ) {
      print FILE '_'.$sample.'('.$$srnas_ref{$seq}{$sample}.')' ;
    }
    print FILE "\n$seq\n" ;
  }
  close FILE ;
  
} #write_nr_seq_file


1;