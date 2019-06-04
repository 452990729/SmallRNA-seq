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


package SrnaTools::Job::Target ;
use base SrnaTools::Job ; 

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
   
  # find targets
  # we don't filter the sRNA file because we need
  # to preserver the IDs in this case but
  # we do uncompress if necessary (standard)
  $self->_run_module(
    'target',{
      transcriptome_file => $self->param('transcriptome'),
      transcriptome_name => $self->get_transcriptome_dname($self->param('transcriptome')),
      infile             => 'srna_file.uncompressed', 
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
  
  # We must have an upload file or some pasted sequences
  if (!$self->param('srna_file') && !$self->param('pasted_srnas') ) {
    $self->_set_param_err('No sRNA file selected for upload and no sRNA sequences pasted ','srna_file','pasted_seq');
  } elsif ($self->param('srna_file') && $self->param('pasted_srnas') ) {
    $self->_set_param_err('Please either select a file to upload or paste sRNA sequences into the textbox.','srna_file','pasted_seq');
  } 
  $self->auto_validate_params('pasted_srnas');
  
} # _do_validate_user_parameters


#######################################################
# Write pasted sequences to file srna_file.raw, if we
# have any, after the job dir is created
#######################################################
sub setup_job_directory{
  my $self = shift;
  my $job_storage_dir_server = shift || SrnaTools::Exception->throw('missing parameter job_dir_eerver') ;
  $self->SUPER::setup_job_directory($job_storage_dir_server);
  
  if ( $self->param('pasted_srnas') ) {
    my $dest_file = $self->job_working_dir.'/data/srna_file.raw';
    
    open (SRNA_FILE, ">", $dest_file) or SrnaTools::Exception::Fileupload->throw(message=>"Could not write pasted sRNA sequences to file",log_msg=>"target_file: $dest_file") ;
    
    print SRNA_FILE $self->param('pasted_srnas') ;
    
    close SRNA_FILE ;
  }
  
}

sub transcriptomes_subfolder{ $_[0]->cfg->{transcriptomes_subfolder} }

# The transcriptome display name is available in the config.
sub get_transcriptome_dname{
  my $self = shift ;
  my $file_name = shift;
  my ($transcr) = grep( $_->{fname} eq $file_name , @{$self->cfg->{transcriptomes}} ) ;
  return $transcr->{dname};
}

# DO NOT DELETE
1 ;