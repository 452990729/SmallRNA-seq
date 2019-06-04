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


package SrnaTools::Job::Filter ;
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
      filter_name => 'filter by sequence properties (low-complexity, size-range)'
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
        same_strand => (defined $self->param('trrna_sense')) ? 1 : 0,
      } 
    ) ;
     $self->_run_module(
      'filter_matches',{
        infile  => 'patman.trrna.out',
        filter_name => 'filter t/rRNA',
      } 
    ) ;
  } 
  
  # filter genome matches
  if (defined $self->param('genome') && $self->param('genome') ne 'none') {
    $self->_run_module(
      'run_patman',{
        ref_seq => $self->param('genome'), 
        infile  => 'nr_seqs.infile',
        outfile => 'patman.genome.out',
        status_msg => 'matching to genome',
      } 
    ) ;
     $self->_run_module(
      'filter_matches',{
        infile  => 'patman.genome.out',
        filter_name => 'filter genome',
        filter_nomatch => (defined $self->param('discard_genome_matches')) ? 0 : 1,
      } 
    ) ;
  }
  
  # build a string describing what we did (for the result file)
  # List the sequences we removed
  my @filters;
  push @filters, 'low complexity';
  if (defined $self->param('trrna')){
    my $filter = 'matching t/rRNAs';
    $filter .= ' (sense strand only)' if (defined $self->param('trrna_sense'));
    push @filters, $filter;
  }
  # the min/max size is empty if '--' was selected (default)
  if ($self->param_or_default('minsize')){
    push @filters, '<'.$self->param_or_default('minsize').'nt';
  }
  if ($self->param_or_default('maxsize')){
    push @filters, '>'.$self->param_or_default('maxsize').'nt';
  }
  if (defined $self->param('genome') && $self->param('genome') ne 'none'){
    my $filter = (defined $self->param('discard_genome_matches')) ? 'matching' : 'not matching';
    $filter .= ' genome '. $self->param('genome');
    push @filters, $filter;
  }
  my $filters_string = join ('; ', @filters);
  
  $self->_run_module(
    'filter_generate_result',{
      infile  => 'nr_seqs.infile',
      make_nr => (defined $self->param('make_nr')) ? 1 : 0,
      filters_string => $filters_string,
    } 
  ) ;
  
 #print Data::Dumper->Dump([$self]);exit;
   

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
  
  # We must have an upload file
  if (!$self->param('srna_file') ) {
    $self->_set_param_err('No FASTA file selected for upload','srna_file');
  }
  # The auto_validate_param method
  # requires a detailed definition of
  # parameters in tools config with
  # def,form,desc etc. The name of the
  # parameter must match the name of the
  # config item. See SrnaTools::Job::auto_validate_param
  # for more details. It is always possible to
  # manually validate everything here as well
  $self->auto_validate_params(
    'minsize',
    'maxsize',
  );
  
  if ($self->param('minsize') && $self->param('maxsize')) {
    if ( $self->param('minsize') > $self->param('maxsize') ) {
      $self->_set_param_err(
        "Maximum sRNA size must be greater than minimum.",
        'maxsize',
        'minsize'
       );
    }
  }
    
  
} # _do_validate_user_parameters


# DO NOT DELETE
1 ;