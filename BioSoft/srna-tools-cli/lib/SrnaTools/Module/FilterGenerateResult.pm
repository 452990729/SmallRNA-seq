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
# It generates the result files from the filtered input
# data set 
#
# Arguments
# --working_dir             Working-directory for this project
# --proj                    Name of job for output
# -i INPUTFILE              (filtered) sRNA input file, assumed to be in data/
# --filters_string          Names of filters used (for overview table)
# --make_nr                 Make output FASTA file non-redundant
#
# We generate the output in the subdir "results" in the working directory. 

package SrnaTools::Module::FilterGenerateResult ;
use base SrnaTools::Module ;
use strict ;
use warnings;

sub run{
my $self = shift ;

#############################
# Parameters and
# declarations
#############################
my $module_name = "filter results writer" ; # for error log
my $working_dir = $self->job->job_working_dir ;
my $job_name = $self->job->job_name;
my $data_dir ;
my %output_opts ;
my $infile_name = $self->param('infile');
my $infile_full_path;
my $filters_string = $self->param('filters_string') || '';
my $make_nr = $self->param('make_nr');


#############################
# Check parameters and files
#############################

eval {
  die "missing working directory parameter - please contact administrator\n" unless $working_dir ;
  die "missing job name parameter - please contact administrator\n" unless $job_name ;
  die "missing input file name\n" unless $infile_name ;
  
  # Default error file is WORKING_DIR/errors
  $data_dir = $working_dir.'/data/' ;
  $infile_full_path = $data_dir.'/'.$infile_name ;
  
  # Check files and directories
  die "could not find/read working directory - please contact administrator\n" unless -d $working_dir ;
  die "data directory not found in working directory - please contact administrator\n" unless -d $data_dir ;  
  
  die "could not find zip executable - please contact administrator\n" unless $self->binary_in_path('zip'); ;
  
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error in parameters: $@\n") ;
}

# Set options for output
$output_opts{job_name} = $job_name ;
$output_opts{filters_string} = $filters_string;
$output_opts{make_nr} = $make_nr;

# update status on server
$self->update_status("generating filter result files") ;


#############################
# Generate output
#############################
eval{
  my $dest_dir = $working_dir.'/results/filter_results' ;
  mkdir $dest_dir or die "could not generate result directory\n" ;
  my $dest_full_path = $dest_dir.'/'.$output_opts{job_name}.'_filtered.fasta' ;
  
  # overview of filtering results
  $self->generate_table(
    $dest_dir,
    \%output_opts,
  ) ;
  
  # process fasta file
  $self->generate_fasta(
    $infile_full_path, 
    $dest_full_path,
    \%output_opts,
  ) ;
  
  # zip results: -q quite, -r recursive, -j exclude full path to dir
  my $cmd = "zip -qrj $dest_dir $dest_dir" ;
  system($cmd) == 0 or die "Failed to compress result files\n" ; 
  system("rm -rf $dest_dir") ; # remove original directory
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while generating result files: $@\n") ;
}
return 1;
} # run


#################### subroutines

#############################
# Generate csv table file
# with overview of filtering steps
#############################
sub generate_table{
  my $self = shift;
  my $dest_dir = shift ;
  my $output_opts_ref = shift ;
  
  my $file_name = $$output_opts_ref{job_name}.'_overview.csv' ;
  open (RESULT , '>', $dest_dir.'/'.$file_name) or die "Could not generate overview file\n" ;
  
  print RESULT "\"Overview of filtering results for $$output_opts_ref{job_name}\"\n" ;
  print RESULT "\"Removed sequences: $$output_opts_ref{filters_string}\"\n" ;
  
  # print table of read counts
  # stages are "input", "filtered: sequence properties", "filtered: t/rRNA"
  # We only have one sample S1 at the moment
  print RESULT "\"Read counts\"\n";
  print RESULT "\"\",\"total\",\"non-redundant\"\n" ;
  foreach my $stage ( @{$self->job->cache->{read_counts}} ) {
    my ($stage_name) = (keys %$stage) ; # there is only one key 
    print RESULT "\"$stage_name\"" ;
    print RESULT ','.$$stage{$stage_name}{'S1'}{total};
    print RESULT ','.$$stage{$stage_name}{'S1'}{nr} ;
    print RESULT "\n" ;
  }
  
  close RESULT ;
} # generate_table

# generate final fasta file
# The filtered file is in nr form.
# if the "make_nr" option for the output
# is selected then simply replace the ID
# lines with consecutive number and ABUNDANCE in
# brackets. Otherwise, print the sequence
# ABUNDANCE times. At this point it should be
# safe to assume that we have ID lines in the
# format 
# >SEQUENCE_S1(ABUNDANCE), so we can get
# both sequence and abundance from the ID
sub generate_fasta{
  my $self = shift;
  my $infile_full_path = shift ;
  my $dest_full_path = shift ;
  my $output_opts_ref = shift ;

  open (IN, '<', $infile_full_path) or die "Could not read sequence file\n" ;
  open (OUT, '>', $dest_full_path) or die "Could not generate sequence file\n" ;

  my $i = 1;
  if ($$output_opts_ref{make_nr}) {
    while (<IN>) {
      if (s/>(\w+?)_[a-z0-9_]+/>$i/) {++$i;}
      print OUT;
    } # while
  } else { # convert to redundant form
    while (<IN>) {
      next unless substr($_,0,1) eq '>' ; # get everything from ID line
      die "could not parse FASTA ID line '$_'" unless />(\w+?)_[a-z0-9_]+\((\d+)\)/i ;
      foreach my $j (1..$2){
        print OUT '>'.$i.'_'.$j.'('.$2.")\n".$1."\n" ;
      }
    ++$i ;
    } # while
  } #if make_nr
  close IN;
  close OUT;

  return ;
}

  
1;