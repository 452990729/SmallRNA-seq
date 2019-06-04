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
# This is the base class for all Modules if the SrnaTools app
# Modules are the parts that do the actual data-analysis work 
# inside the App, they are run by subclasses of Job.pm 
#
#####################################################################

package SrnaTools::Module ;
use strict;
use warnings;
use SrnaTools::Exceptions ;
use Exception::Class::TryCatch;

#######################################################
# A new instance needs a reference to the Job object 
# that is calling it.
# It may take an additional hash ref of arguments
#######################################################
sub new{
  my ($classname, $job_obj, $params ) = @_ ;
  SrnaTools::Exception->throw("missing argument: job_object") unless $job_obj;
  SrnaTools::Exception->throw("argument job_object is not of SrnaTools::Job class.") unless $job_obj->isa("SrnaTools::Job");
  my $self = {
    _job    => $job_obj,
    _params => $params,
  } ;
  bless $self, $classname;
  
  # In some environments, CPAN modules might be in 
  # non-standard directories. 
  # Affected modules must be included by a "require"
  # statement in the 'run' method of the module. Here
  # we add to @INC if additional paths defined in config
  # but we can't use "use lib" because this happens at
  # run time, not compile time
  my $module_dirs;
  eval{
    $module_dirs = $self->job->app->path_to('module_dir',$self->job->execution_env) ;
  };
  if ($module_dirs && ref $module_dirs eq 'ARRAY'){ 
    @INC = (@$module_dirs, @INC);
  }
  return $self ;
}
#######################################################
# accessors
#######################################################
sub job{ $_[0]->{_job} }

# module parameters (passed from Job)
# $self->params - hash ref of all params
# $self->param('name') - value of named parameter
sub params{ $_[0]->{_params} }
sub param{
  my ($self, $param_name) = @_ ;
  SrnaTools::Exception->throw("Missing argument in sub: no parameter name given") unless $param_name;
  return $self->params ? $self->params->{$param_name} : undef ;
}

#######################################################
# Delegate update_status calls to the instance of the
# App that is running this module
#######################################################
sub update_status{
  my $self = shift ;
  my $status_msg = shift ;
  $self->job->app->update_job_status($self->job, $status_msg);
}


#######################################################
# check that a tird party binary file is in the PATH
#######################################################
sub binary_in_path{
  my ($self, $binary) = @_ ;
  
  my $out = `which $binary`;
  if (!$out || $out=~/which: no/){
    return 0 ;
  } else {
    return 1;
  }
}

# DO NOT DELETE
1;