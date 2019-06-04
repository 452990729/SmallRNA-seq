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


package SrnaTools::Job::Siloco ;
use base SrnaTools::Job ; # requires Tool and isa Tool

use strict ;
use warnings ;

#######################################################
# The constructor is only implemented 
# in the base class
#######################################################


#######################################################
# Accessors
#######################################################

#######################################################
# accessors for tool config items
#######################################################


#######################################################
# Some additional tool-specific initialisation
# of new objects.
# Here we need to get the name of the mirbase file
# which depends on the setting of 'overhangs'
#######################################################
sub _init{
  my $self = shift ;
  
  # First do whatever the base class wants
  $self->SUPER::_init(@_) ;
  
}


#######################################################
# Run the job.
# The code that carries out the actual data preparation
# and analysis is organised into Modules (see 
# lib/SrnaTools::Module).
#
# Modules have full access to the Job and App instances
# that run them, so they can access config data directly.
# Parameters could also be accessed directly but they
# are better handled here, so we explicityl pass them
# in using the _run_module method:
# $self->_run_module('name_of_module', {parameters})
#
# The module name can be given in lower case with
# underscores or in camel-case, so 
# 'my_module' and 'MyModule' would both load a package
# lib/SrnaTools/Module/MyModule.pm
#######################################################
sub run{
  my $self = shift ;
  
  $self->SUPER::run(@_);
   #print Data::Dumper->Dump([$self]);exit;
   
  $self->_run_module(
    'filter_by_seq_properties',{
      filter_low_comp => 1,
      min_length => $self->param_or_default('minsize'),
      max_length => $self->param_or_default('maxsize'),
      outfile => 'nr_seqs.infile',
    } 
  ) ;

  # filter t/rRNA matches
  if (defined $self->param('trrna') ) {
    $self->_run_module(
      'run_patman',{
        ref_seq => $self->trrna_file, 
        infile  => 'nr_seqs.infile',
        outfile => 'patman.trrna.out',
        status_msg => 'matching to t and rRNAs',
      } 
    ) ;
     $self->_run_module(
      'filter_matches',{
        infile  => 'patman.trrna.out',
        filter_name => 'filter t/rRNA',
      } 
    ) ;
  } 
  
  # get genome matches
  $self->_run_module(
    'run_patman',{
      ref_seq => $self->param('genome'), 
      infile  => 'nr_seqs.infile',
      outfile => 'patman.out',
      status_msg => 'matching to genome',
      ensure_result => 1, # no matches are treated as error
      remove_infiles => 1, # delete raw seq ASAP
    } 
  ) ;
  
  # genome browser links
  my @gb_links ;
  push @gb_links, 'tair' if (defined $self->param('tair_links'));
  push @gb_links, 'asrp' if (defined $self->param('asrp_links'));
  
  # association between file name parameter and
  # sample name for job cache
  $self->cache->{file_to_sample_param} = {
    srna_file1 => $self->param('sample_name1'),
    srna_file2 => $self->param('sample_name2'),
  };

  # run the locus comparison
  $self->_run_module(
    'siloco',{
      min_hits          => $self->param_or_default('min_hits'),
      max_gap           => $self->param_or_default('max_gap'),
      pseudocount       => $self->param_or_default('pseudocount'),
      single_sample     => $self->param_or_default('num_samples') eq '1' ? 1 : 0,
      show_unique       => (defined $self->param('uniq')) ? 1 : 0,
      sort_method       => 'pos', # can also be 'rank' - make optional
      patman_out_file   => 'patman.out',
      ratio_rank_weight => 0.5,
      avg_rank_weight   => 0.5,
      gb_links          => \@gb_links,
    } 
  ) ;
  

} # run


#######################################################
# Validate input parameters and populate error
# attributes. First call the method in Tool parent
# class for some generic validation then do some stuff
# specific for this tool.
#######################################################
sub _do_validate_user_parameters{
  my $self = shift ;
  my $errors = [] ;
  my $cfg = $self->config ;
  
  # Do some generic validation that is
  # required by all tools and implemented
  # in the parent class
  $self->SUPER::_do_validate_user_parameters ;
  
  ####### This is tool specific validation ############
  
  $self->auto_validate_params(
    'sample_name1',
    'sample_name2',
    'num_samples',
    'min_hits',
    'max_gap',
  );
  
  # Depending on the mode we need either one or two data files
  foreach (1..$self->param_or_default('num_samples')) {
    if (!$self->param('srna_file'.$_) ) {
      my $add_txt ='';
      if ($_ == 2) {
        $add_txt = "(or select single-dataset mode)"
      }
      $self->_set_param_err("Upload a FASTA file in file field $_ $add_txt",'srna_file'.$_);
    }
    
    if (!$self->param('sample_name'.$_) ) {
       $self->_set_param_err("Name of sample $_ missing",'sample_name'.$_);
    }
  }
  
  if ( $self->param_or_default('num_samples') == 1) {
    if ($self->param('srna_file2') ) {
      $self->_set_param_err("Single sample mode selected but 2 files chosen for upload.",'srna_file2','num_samples');
    }
  }
  
  # Can't have two identical sample names
  if ($self->param('sample_name1') eq $self->param('sample_name2')) {
       $self->_set_param_err("Names of samples can not be identical",'sample_name1','sample_name2');
  }
  
  # Must have a genome
  my $genome = $self->param('genome') ;
  if (! $genome or $genome=~/please select/ or $genome eq 'none') {
    $self->_set_param_err('No genome selected','genome');
  }
  
  
} # _do_validate_user_parameters



# DO NOT DELETE
1 ;