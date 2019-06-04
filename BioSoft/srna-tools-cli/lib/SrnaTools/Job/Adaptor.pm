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


package SrnaTools::Job::Adaptor ;
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
    'adaptor',{
      infile             => 'srna_file.uncompressed',
      adaptor_sequence_3 => $self->param('adaptor_sequence_3'),
      adaptor_sequence_5 => $self->param('adaptor_sequence_5'),
      allow_revcomp      => (defined $self->param('allow_rev_comp')) ? 1 : 0,
      minsize            => $self->param_or_default('minsize'),
      maxsize            => $self->param_or_default('maxsize'),
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
    'adaptor_sequence_3',
    'adaptor_sequence_5',
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