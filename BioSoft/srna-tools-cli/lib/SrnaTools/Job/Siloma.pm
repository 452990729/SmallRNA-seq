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


package SrnaTools::Job::Siloma ;
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
  
  # this module either has to fetch from genome or write
  # the pasted seq into a file. 
  $self->_run_module(
    'get_backbone_seq',{
      outfile       => 'backbone_seq.fasta',
      refseq        => $self->param('genome'),
      region_chrom  => $self->param('region_chrom'),
      region_start  => $self->param('region_start'),
      region_end    => $self->param('region_end'),
      genome_name   => $self->get_genome_dname,
      pasted_seq_ref => \$self->param('pasted_seq'),
    } 
  ) ;

  # get matches
  $self->_run_module(
    'run_patman',{
      ref_seq => 'backbone_seq.fasta',
      ref_seq_dir => $self->job_working_dir.'/data/',
      infile  => 'nr_seqs.infile',
      outfile => 'patman.out',
      status_msg => 'matching to reference sequence',
      remove_infiles => 1, # delete raw seq ASAP
    } 
  ) ;
  
  # generate result files
  $self->_run_module(
    'siloma',{
      plot_nr       => (defined $self->param('plot_nr')) ? 1 : 0,
      plot_labels   => (defined $self->param('plot_labels')) ? 1 : 0,
      pat_out_file  => 'patman.out',
      backbone_file => 'backbone_seq.fasta', 
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
  
  # We must have an upload file
  if (!$self->param('srna_file') ) {
    $self->_set_param_err('No sRNA file selected for upload','srna_file');
  }
  
  $self->auto_validate_params(
    'minsize',
    'maxsize',
  );
  
  my $genome = $self->param('genome') ;
  my $pasted_seq = $self->param('pasted_seq') ;
  
  if ($genome && $genome ne 'none' && $genome!~/please select/ ){
    $self->auto_validate_params('region_chrom','region_start','region_end');
    if (! $self->param('region_chrom') or
        ! $self->param('region_start') or
        ! $self->param('region_end') ) {
      $self->_set_param_err(
        'Please select chromosome, start and end position of the genomic region','region_chrom','region_start','region_end'
      );
    } else {
      if ( $self->param('region_start') >= $self->param('region_end') ){
        $self->_set_param_err(
        'Start position must be smaller than end','region_start','region_end'
        );
      }
    }
  } elsif ( $pasted_seq ) {
    $self->auto_validate_params('pasted_seq') ;
  } else {
    $self->_set_param_err(
      'No reference sequence given: either select a genome region or paste a sequence into the textbox','genome','region_chrom','region_start','region_end', 'pasted_seq'
    );
  }
  
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