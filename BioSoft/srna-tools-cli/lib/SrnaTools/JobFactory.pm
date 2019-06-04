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
# This is an abstract factory class that instantiates Job objects
# either by receiving a tool name and parameters/config or from an
# existing job directory, where parameters, including tool name, 
# can be obtained from config files.
#
#####################################################################

package SrnaTools::JobFactory ;
use strict ;
use warnings ;

# A new tool object needs a config data structure
# and the parameters for running th tool (e.g. from 
# CGI object or command line)
# The parameters are optional because we also need the 
# object to set up input forms before any parameters
# can be passed to the tool
sub create_job {
  my ($classname, %args) = @_ ;
  my $tool_name = $args{tool_name} ;
  if (!$tool_name) {
    SrnaTools::Exception->throw("missing tool_name argument in subroutine") ;
  }
  
  my $class = "SrnaTools::Job::"._camelcase($tool_name);
  
  # untaint class name (important for web version)
  if ($class!~/^([a-z:]+)$/i){ 
    SrnaTools::Exception->throw(message=>"invalid job class name",log_msg=>"'$class' did not pass regex");
  } else {
    $class = $1 ;
  }
  
  # Now load the module
  eval "require $class" ;
  if ($@) {
    if ($@=~/Can't locate/i) {
      SrnaTools::Exception->throw("Could not find module $class") ;
    } else {
      SrnaTools::Exception->throw(message=>"Module $class has errors",log_msg=>"$@") ;
    }
  }
  
  return $class->new( 
    app             => $args{app},
    tool_name       => $tool_name, 
    tool_conf_file  => $args{tool_conf_file}, 
    params          => $args{params}, 
    app_conf_file   => $args{app_conf_file},
    tool_conf_file  => $args{tool_conf_file},
    params_file     => $args{params_file},
    job_working_dir => $args{job_working_dir}, 
  )
} # new


#######################################################
# (Re-)Create a job from the job config file(s).
# This is used to run a job remotely or through batch 
# queue and to manage remote jobs on the server.
# Arguments:
# job_dir: required path to the working dir of this 
#          job, which must contain a config sub-dir
# optional:
# skip_tool_conf, skip_app_conf, skip_params
# set these to 1 to skip passing file names to the
# Job. This is useful in the WebApp when we just need 
# to check the status of the job and reading these files
# would be a waste of time.
#######################################################
sub create_job_from_config{
  my ($classname, %args) = @_ ;
  my $job_dir = $args{job_dir} ;
  if (!$job_dir) {
    SrnaTools::Exception->throw("missing job_dir argument in subroutine") ;
  }
  
  my $cfg_dir = $job_dir.'/config';
  if (! -r $cfg_dir) {
    SrnaTools::Exception::FileAccess->throw(message=> "Could not access job config directory", log_msg=>"Trying to recreate job from the config directory '$cfg_dir' but it is not present or not readable");
  }
  
  # Fetch tool name from file in 'config'
  my $tool_name_file = $cfg_dir.'/tool_name' ;
  if (! -r $tool_name_file) {
    SrnaTools::Exception::FileAccess->throw(message=> "Could not access tool name file in job config directory", log_msg=>"trying to access: $tool_name_file");
  }
  
  # Get name of tool from file and untaint
  # it (will be checked further in Job class)
  ($tool_name_file) = ($tool_name_file=~/^(.+)$/) ;
  my $out = `cat $tool_name_file 2>/dev/null` ;
  my ($tool_name) = ($out=~/^(.+)$/) ;
  if (! $tool_name) {
     SrnaTools::Exception->throw(message=>"Could not read tool name from config",log_msg=>"no output from 'cat $tool_name_file '") ;
  }
  my $job_args = { tool_name => $tool_name } ;
  my $tool_conf_file = $args{skip_tool_conf}  ? '' : $cfg_dir.'/tool.conf';
  my $app_conf_file = $args{skip_app_conf}  ? '' : $cfg_dir.'/application.conf';
  my $params_file = $args{skip_params} ? '' : $cfg_dir.'/params';
   
  return SrnaTools::JobFactory->create_job( 
    app             => $args{app},
    tool_name       => $tool_name,
    app_conf_file   => $app_conf_file,
    tool_conf_file  => $tool_conf_file,
    params_file     => $params_file,
    job_working_dir => $job_dir,
  ) ;
}



sub _camelcase{
  my $name = shift ;
  $name =~s/(\b|_)([a-z])/\u$2/ig ;
  return $name ;
}


# DO NOT DELETE
1 ;