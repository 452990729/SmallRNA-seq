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
# This is the parent class for all sRNA tool jobs
#
#####################################################################

package SrnaTools::Job ;
use strict ;
use warnings ;
use File::Temp qw( tempfile tempdir );
use File::Basename;
use Data::Dumper ;
use Template ;
use SrnaTools::Exceptions;

#######################################################
# A new Job object needs a tool config 
# but the (user-submitted) run-time parameters are not 
# obligatory because we also create Job objects to 
# render the input form before any parameters can
# be submitted. 
#######################################################
sub new {
  my $classname = shift ;
  my $self = {};
  bless $self, $classname;
  $self->_init(@_) ; # run additional initialisation
  return $self ;
} # new


#######################################################
# Run the job. 
# We get a connection to the App instance that called
# this by passing in an App object. 
# The app config is updated to use the one stored with
# this job (the App itself does that).
# TODO: turn all the back-end scripts into
# SrnaTools::Module classes so that they can do proper
# error handling and be integrated into Jobs by
# calling their run method from here. This would also
# make the status updates a lot easier.
# Problem: tricky to 'require' the correct package
# files if we want to have a development and a production
# site on the same machine. 
#######################################################
sub run{
  my $self = shift ;

  #if ($self->param('do_uncompression') ){
    $self->_run_module('uncompress_infiles',{} )
  #}

  # Implement the actual pipeline in tool-specific sub classes

} # run


#######################################################
# Run a module/script from the tools backend. These
# are the modules that perform the actual job.
# Usage:
# $job->_run_module('NAME_OF_MODULE',{ HASH OF OPTIONAL ARGS })
# The name of the module is either the real package name
# (in camel-case) or the name of the script without the
# .pl suffix, e.g.:
# MyModule or my_module
# both look for package SrnaTools::Module::MyModule
# A full abs path to the lib dir for the job's execution
# environment is used to load the module.
#######################################################
sub _run_module{
  my $self = shift ;
  my $module = shift ;
  my $module_params = shift ;
  my $machine =  $self->execution_env ;
  
  SrnaTools::Exception->throw(message=>"No module name given") unless $module ;
  if ($module!~/[A-Z]/){ # need to camel-case module name
    $module = $self->_camelcase($module) ;
  }
  
  # Load the package with a full abs path to 
  # the lib directory for this job's environment
  my $class = "SrnaTools::Module::$module";
  my $class_path = $class ;
  $class_path =~s/::/\//g;
  my $lib_dir = $self->app->path_to('lib_dir',$self->execution_env);
  my $full_class_path = $lib_dir."$class_path.pm";

  # Now load the module
  eval "require '$full_class_path'" ;
  if ($@) {
    if ($@=~/Can't locate/i) {
      SrnaTools::Exception::FileAccess->throw(message=>"Could not find or load package $class or the package has compilation errors.",log_msg=>"$@") ;
    } else {
      SrnaTools::Exception::FileAccess->throw(message=>"Package $class has errors",log_msg=>"$@") ;
    }
  }
  
  if ($module_params && ref $module_params ne 'HASH') {
    SrnaTools::Exception::FileAccess->throw(message=>"module parameters not a hash reference",log_msg=>"usage of _run_module: $self->_run_module('name_of_module',{param1=>value1,param2=>value2,...})") ;
  }

  my $module_obj = $class->new( $self, $module_params ) ;
#   if (! $module_obj->run ) {
#     SrnaTools::Exception::ModuleExecution->throw(message=>"Unknown exception while running module $module", log_msg=>"method 'run' called on $module returned false. The script may have terminated prematurely without throwing a proper SrnaTools::Exception. Also check that the 'run' method returns true as the last statemet in package $module");
#   }
  return $module_obj->run ;
  
} # _run_module


#######################################################
# Initialisation
#######################################################
sub _init{
  my $self = shift ;
  my %args = @_ ;
  
  my $tool_name = $args{tool_name} ;
  SrnaTools::Exception->throw("missing argument: tool name") unless $tool_name;
  if ($tool_name!~/^[\w]+$/) {
    SrnaTools::Exception->throw(message=>"invalid tool name",log_msg=>"Tool names as defined in tool config should only contain letters, number and underscores. The malformed tool name was '$tool_name'. Change tool name in tools config file.") ;
  }
  $self->{_tool_name} = $tool_name ;
  
  # A job can have a connection back to an
  # instance of the app that called it
  if ( $args{app} ) {
    if ( ! $args{app}->isa('SrnaTools::Application') ) {
      SrnaTools::Exception->throw(message=>"invalid instance of Application object",log_msg=>"value of app was: ". ref $args{app}. " and not an instance of SrnaTools::Application");
    }
    $self->{_app} = $args{app};
  }
  
  $self->{_environment} = $args{environment} if $args{environment} ;
  
  # Try to get job parameters (from user)
  # either from passed-in argument or from
  # a file.
  if ( $args{params} ){
    $self->{_params} = $args{params} ;
  } elsif ( $args{params_file} ) {
    $self->{_params} = $self->_parse_config_file( $args{params_file} ) ; 
  }
  
  # Parse tool and app config from 
  # files if given
  if ( $args{tool_conf_file} ) {
    my $cfg_ref = $self->_parse_config_file( $args{tool_conf_file} ) ;
    $cfg_ref = $self->_merge_tool_conf($cfg_ref) ;
    $self->{_cfg} = $cfg_ref ;
  }
  
  # NOTE should not need to load app config here. Just store
  # path to it so the job can reconfigure the running app 
  if ( $args{app_conf_file} ) {
    #$self->{_application_config} = $self->_parse_config_file( $args{app_conf_file} ) ;
    $self->{_application_config_file} = $args{app_conf_file} ;
  }
  
  # Set the job name from user input or
  # use the MyJob default.
  my $job_name = $self->param('job_name') || 'MyJob' ;
  SrnaTools::Exception->throw("Can not create job without a name") unless $job_name ;
  $job_name=~s/\s/_/g ; # job name is used as file name, so can't have spaces
  $self->job_name($job_name) ;
  
  # A new job generates its own working directory
  # but if we "re-constitute" an already prepared job
  # we should take its path as the working dir
  # The JobFactory handles this automatically
  $self->job_working_dir( $args{job_working_dir} ) if $args{job_working_dir} ;
  
  # Jobs have a cache that can hold data between modules
  # such as the number of filtered sequences
  $self->{_cache} = {} ;
  
} # _init


#######################################################
# Read config or parameter files using Config::Auto
#######################################################
sub _parse_config_file{
  my ($self, $cfg_file) = @_ ;
  if (! $cfg_file) {
    SrnaTools::Exception->throw("missing argument: cfg_file") ;
  }
  if (! -r $cfg_file) {
    SrnaTools::Exception::FileAccess->throw(message=> "Could not access tool config directory", log_msg=>"trying to access: $cfg_file");
  }
  ($cfg_file) = ($cfg_file=~/(.+)/) ; # untaint
  my $cfg_ref ;
  eval {
    $cfg_ref = Config::Auto::parse( $cfg_file, format => "perl" );
  } ;
  if ($@) {
    SrnaTools::Exception::FileAccess->throw(message=> "Could not read tool config file", log_msg=>"trying to access: $cfg_file. Got: $@") ;
  }
  
  return $cfg_ref ;
  
}

#######################################################
# Merge tool config file sections in such a way that
# tool-specific config overwrites the defaults given
# in the 'all' section of the tools config file.
#######################################################
sub _merge_tool_conf{
  my ($self, $cfg_ref) = @_ ;
  my $tool_name = $self->tool_name ;
  
  if (! $tool_name) {
    SrnaTools::Exception->throw(message=>"could not retrieve tool name for merging tool config data",log_msg=>"_merge_tool_conf method was called and it assumes that the tool has a name already but no tool name could be found in instance var _tool_name. Something went wrong whe ninitialising the Job object") ;
  }
  if (! $cfg_ref) {
    SrnaTools::Exception->throw("missing argument: cfg_ref") ;
  }
  
  # If we are reading a copy of the tools config
  # then it has already been merged and we don't 
  # have these sections
  return $cfg_ref if ! defined $cfg_ref->{$tool_name} ;
  
  # merge section for this tool on
  # top of 'all' by evaluating both as
  # a single list an putting back into 
  # a hash. 
  $cfg_ref = { ( %{$cfg_ref->{all} }, %{$cfg_ref->{$tool_name} } ) } ;
  
}


#######################################################
# Accessor methods
#######################################################
sub tool_name{ $_[0]->{_tool_name} }
sub config{ $_[0]->{_cfg} } # access the tool config
sub cfg{ $_[0]->{_cfg} } # alias for config
sub cache{  $_[0]->{_cache} }
sub is_web { $_[0]->{_environment} eq 'web' ? 1 : 0 }
sub app_config_file{ $_[0]->{_application_config_file} }
sub app{
  $_[0]->{_app} || SrnaTools::Exception->throw("No App object stored");
}

# The job name is given by the user, The ID is a unique 
# tag that is generated by letting tmpdir generate a 
# directoy name by concatanating the job name and a 
# random string. These are just accessors to both.
# See setup_job_directory for details.
# Both job_name and id are also stored in the params of a job
# so try to get them from there if not stored already in instance var.
sub job_name{
  my ($self,$name) = @_;
  $self->{_job_name} = $name if $name ;
  if (! defined $self->{_job_name}) {
    $self->{_job_name} = $self->param('job_name') ;
  }
  return $self->{_job_name};
}
sub job_id{ 
  my ($self,$id) = @_;
  $self->{_job_id} = $id if $id ;
   if (! defined $self->{_job_id}) {
    $self->{_job_id} = $self->param('job_id') ;
  }
  return $self->{_job_id};
}

# This is the complete abs paths to the working 
# directoy of this job on the local machine
# It is either set by setup_job_directory
# upon creation of a job or passed in as
# arg if we re-create an existing job from files
sub job_working_dir{
 my ($self,$dir) = @_;
  $self->{_job_working_dir} = $dir if $dir ;
  return $self->{_job_working_dir};
}


#######################################################
# Access tool parameters.
# List of all user-parameters:
# $tool->params
# Value of specific parameter (retrieve or set):
# $tool->param('param-name'[, new_value])
#######################################################
sub params{ 
  my $self = shift ;
  shift && SrnaTools::Exception->throw(message=>"wrong number of arguments passed to subroutine",log_msg=>"method 'params' takes no arguments - use 'param(NAME) instead") ;
  $self->{_params} ;
}
sub param{
  my ($self, $param_name, $new_value) = @_ ;
  SrnaTools::Exception->throw("Missing argument in sub: no parameter name given") unless $param_name;
  return unless $self->params ;
  $self->params->{$param_name} = $new_value if $new_value ;
  return  $self->params->{$param_name} ;
}
# Like param but attempts to retrieve the default
# if given in tools config like this:
# param_name => {
#   desc=>"A parameter",
#   form=>'i', # an integer
#   def=> 10, # the default
# }
sub param_or_default{
  my ($self, $param_name) = @_ ;
  SrnaTools::Exception->throw("Missing argument in sub: no parameter name given") unless $param_name;
  return unless $self->params ;
  if (defined $self->params->{$param_name} ) {
    return $self->params->{$param_name};
  } elsif (defined $self->cfg->{$param_name}{def} ){
    return $self->cfg->{$param_name}{def};
  } else {
    return undef ;
  }
}

#######################################################
# Accessors for common config items 
#######################################################
# direct access
sub genomes{ $_[0]->cfg->{genomes} }
sub trrna_file{ $_[0]->cfg->{trrna_file} }
sub max_number_uniq_seqs{ $_[0]->cfg->{max_number_uniq_seqs} }
sub display_name{ $_[0]->cfg->{display_name} || $_[0]->tool_name || ''}
sub description{ $_[0]->cfg->{description} }
sub requires_email{$_[0]->cfg->{requires_email} }
sub min_seq_size{ $_[0]->cfg->{min_seq_size} }
sub queue{ $_[0]->cfg->{queue} }
sub result_page_partial{ $_[0]->cfg->{result_page_partial} }

sub mode{ 
  my ($self,$value) = @_;
  $self->cfg->{mode}=$value if $value;
  return $self->cfg->{mode};
}

# return an array of parameter field names
# that should hold data files
sub data_file_params{
  my $self = shift ;
  my $fields = $self->cfg->{data_file_params} ;
  return unless $fields ;
  $fields = [ $fields ] unless ref $fields eq 'ARRAY';
  return $fields ;
}

# Files which are rendered into templates
# These contain tool logos and names.
# Use 'TOOLNAME_form.tt' and
# TOOLNAME_head.tt' as defaults
sub form_partial{
  my $self = shift ;
  if ($self->cfg->{form_partial} ){
    return $self->cfg->{form_partial} ;
  } else {
    return $self->tool_name.'_form.tt' ;
  }
}
sub head_partial{
  my $self = shift ;
  if ($self->cfg->{head_partial} ){
    return $self->cfg->{head_partial} ;
  } else {
    return $self->tool_name.'_head.tt' ;
  }
}

# derived information

# Some information about the run mode, which can be
# instant, queue_local or queue_remote. Provide some shortcuts
#sub run_remotely{ $_[0]->cfg->{mode}=~/remote/ ? 1 : 0 }
sub queue_mode{ $_[0]->cfg->{mode}=~/queue/ ? 1 : 0 }
sub instant_mode{ $_[0]->mode eq 'instant' ? 1 : 0 }
sub queue_remote_mode{ $_[0]->mode eq 'queue_remote' ? 1 : 0 }
sub queue_local_mode{ $_[0]->mode eq 'queue_local' ? 1 : 0 }
sub remote_mode{ $_[0]->mode=~/remote/ ? 1 : 0 }

# The execution environment is used to select
# paths etc. according to where the job will
# be executed. server and local_queue refer
# to the same machine but they may have different
# lib directories etc.
sub execution_env{ 
  $_[0]->queue_remote_mode ? 'cluster' : 
  $_[0]->queue_local_mode ? 'local_queue' :
  'server' 
}

#######################################################
# Error handling for user input, i.e. errors that should
# not lead to termination of the program.
# $job->param_err_msgs
# returns an array of validation error messages
# $job->param_err_names
# returns an array of input field names with errors if
# possible. These can be used in a WebApp to highlight
# the fields that need to be corrected.
#######################################################
sub param_err_msgs{ return $_[0]->{_param_err_msgs} }
sub param_err_names{ return $_[0]->{_param_err_names} }

# The genome display name is available in the config.
# If we are in CLI App interactive mode then the genome_dname
# paramter is set by the App and we return that
sub get_genome_dname{
  my $self = shift ;
  return $self->param('genome_dname') if $self->param('genome_dname');
  my $genome_file = $self->param('genome');
  if ( !$genome_file || $genome_file eq 'none') {
    return 'none' ;
  }
  my ($gn) = grep( $_->{fname} eq $genome_file , @{$self->cfg->{genomes}} ) ;
  return $gn->{dname};
}


#######################################################
# Validate user input (parameters for the job). Return 1 if
# all ok, and 0 if there are errors. Also set a flag
# once validation is done. This method calls the private
# method _do_validate_user_parameters which populates error
# attributes of the object that can then be retrieved 
# with get_param_err_msgs and get_param_err_names.
#######################################################
sub validate_user_parameters{
  my $self = shift ;
  
  # do the actual error checking unless we
  # have done it already
  unless ($self->{_params_validated} ) {
    $self->_do_validate_user_parameters ;
    $self->{_params_validated} = 1 ;
  }
  
  # validation was successful if we have no errors now
  return $self->has_param_errs ? 0 : 1 ;
  
} # validate_user_parameters


#######################################################
# ask if we have parameter validation errors
#######################################################
sub has_param_errs{
  my $self = shift ;
  return $self->{_param_err_msgs} ? 1 : 0 ;
}


#######################################################
# Some generic validation of user parameters
# that is common to all tools. The tool-specific
# stuff is implemented in the child class
#######################################################
sub _do_validate_user_parameters{
  my $self = shift ;
  my $errors = [] ;
  
  # If email is a required parameter, make sure that
  # we have one and that it passes some basic validation
  if ($self->requires_email && !$self->params->{email}) {
    $self->_set_param_err("An email address is required for this tool",'email') ;
  }
  if ($self->params->{email} && ! $self->_email_is_valid($self->params->{email}) ) {
    $self->_set_param_err("Email address is invalid. Some unusual addresses may also cause this error. Please enter an alternative address in that case or contact us.",'email') ;
  }
  
  if ($self->app->environment eq 'cmd' && (!$self->param('out') && !$self->param('outdir') ) ){
    $self->_set_param_err("Parameter --out|outdir (target for result files) is required in command line mode", 'out');
  }
  
  # The rest ot the validation should be done
  # in the tool-specific subclasses of Job
  
} # _do_validate_user_parameters


#######################################################
# Validate parameter format. Format can be given
# in tools config like this:
# my_param => {format => X [, default => Y,min => Z...}
# where X should be a string containing one of:
# 'i'   positive integer
# 'f'   positive float
# 's'   string
# '/P/' regex, where P is the pattern
# 'enum' array of possible (string) values must be
#        given as 'vals' or 'values'
#
# further options:
# blank => 'string' to define a value that is trated
# as undefined, e.g. '--' for web forms
# use:
# not_null => 1 to enforce a value for this parameter
# use:
# form_desc =>'SOMETHING'
# To provide a description of the required format that
# replaces the standard text. This is useful for 
# regular expressions (if not given, the regex is shown)
# use:
# max_len => x
# to define a maximum length for regex or string 
#
# Errors will be registered as parameter err.
# Will throw an error if the item does not exist 
# or no format has been defined.
# Can also use a list of params like this:
# $self->auto_validate_params(
#    'minsize',
#    'maxsize',
#    'min_abundance'
# );
#######################################################
sub auto_validate_params{
  my $self = shift ;
  foreach (@_) {
    $self->auto_validate_param($_);
  }
}
sub auto_validate_param{
  my $self = shift ;
  my $param = shift ;
  my $format ;
  
  SrnaTools::Exception->throw("No parameter name given") unless $param ;
  
  my $value = $self->param($param) ;
  
  # currently we are not allowing list
  # parameters. Treat this an an excpetion
  if (ref $value eq 'ARRAY' or ref $value eq 'HASH'){
    SrnaTools::Exception->throw("Parameter $param contained a ".ref $value." - list parameters are not allowed and there is probably something wrong with the input form.") ; 
  }
  
  if ( ! $self->cfg->{$param} ) {
    SrnaTools::Exception->throw("can't validate parameter format: no config item with name '$param' exists in tools config") ;
  } elsif ( ! $self->cfg->{$param}{form} )  {
    SrnaTools::Exception->throw("can't validate parameter format: no format ('form') has been defined for config item '$param' in tools config") ;
  } else {
    $format = $self->cfg->{$param}{form} ;
  } 
  
  my $description = $self->cfg->{$param}{desc} || $param ;
  
  # to allow e.g. '--' from a webform to be ignored
  # and treated as a blank value
  # if (defined $self->cfg->{$param}{blank} ) {
  #  $value = undef if $value eq $self->cfg->{$param}{blank} ;
  # }
  
  # If no parameter of this name has been given 
  # return true unless not_null is set 
  # Then return
  if (! defined $value or $value=~/^\s*$/){
    if ($self->cfg->{$param}{not_null} ) {
      $self->_set_param_err(
        "$description: a value is required", $param
      ) ;
    } 
    return; 
  }
  
  # optional description of format from 
  # config. Useful for regex format
  my $format_desc = $self->cfg->{$param}{form_desc} || '' ;
  
  # We have a value - check it's format
  if ($format eq 'i') { # it's an integer
    if ($value!~/^-?\d+$/){
      $format_desc ||= 'must be an integer';
      $self->_set_param_err(
        "$description $format_desc", $param
      ) ;
    } else {
      my $min = $self->cfg->{$param}{min} ;
      my $max = $self->cfg->{$param}{max} ;
      if (defined $min and $value < $min){
        $self->_set_param_err(
          "$description must be >= $min", $param
        ) ;
      } elsif (defined $max and $value > $max){
        $self->_set_param_err(
          "$description must be <= $max", $param
        ) ;
      }
    }
    
  } elsif ($format eq 'f' ) { # a float
    if ($value!~/^-?\d+\.?\d*$/) {
      $self->_set_param_err(
        "$description must be a floating point number", $param
      ) ;
      
    } else {
      my $min = $self->cfg->{$param}{min} ;
      my $max = $self->cfg->{$param}{max} ;
      if (defined $min and $value < $min){
        $self->_set_param_err(
          "$description must be >= $min", $param
        ) ;
      } elsif (defined $max and $value > $max){
        $self->_set_param_err(
          "$description must be <= $max", $param
        ) ;
      }
    }
    
  } elsif ($format eq 's' ) { #string: test only for max length
    my $max_len = $self->cfg->{$param}{max_len} ;
    if (defined $max_len && length($value) > $max_len ){
      $self->_set_param_err(
        "$description must be less than $max_len characters long.", $param
      ) ;
    }
    
  } elsif ($format=~/^\/(.*)\/(i?)$/) { # it's a regex
    my $regex= $1;
    my $mod = $2;
    unless ( $value=~/$regex/ or ($value=~/$regex/i && $mod eq 'i') ) {
      $format_desc ||= 'is not in the correct format (please see manual)';
      $self->_set_param_err(
        "$description $format_desc", $param
      ) ;
    }
    
  } elsif ($format eq 'enum') { # we have an array of allowed values
    my $values = $self->cfg->{$param}{vals} || $self->cfg->{$param}{values};
    if (! $values or ref $values ne 'ARRAY') {
      SrnaTools::Exception->throw("can't validate parameter: no array of string values is given for enum parameter $param.") ;
    }
    if (! grep($_ eq $value, @$values)){
       $self->_set_param_err(
        "$description: value must be one of: ".join(',',@$values), 
        $param
      ) ;
    }
  } else {
    SrnaTools::Exception->throw("can't validate parameter format: don't recognize the format: '$format'") ;
  }
  
}

#######################################################
# Register a parameter error. These are errors that
# arise from validation of parameters against the tool's
# expectations and can be used to generate a list of
# incorrect parameter values on the web form or the
# command line. Each error has to consist at least
# of a string (the error message). Any additional
# arguments are treated as a list of names of the 
# associated parameters that need to be highlighted.
# Example: parameters min = NUM1, max = NUM2, if
# NUM1 > NUM2 we want to throw an error that highlights
# both parameters (e.g. in a web form) and creates one
# error message:
# $self->_set_param_err("Minimum must be less than maximum",
# 'min','max')
# The parameter names are also associated with messages
#######################################################
sub _set_param_err{
  my $self = shift ;
  my $message = shift ;
  my @param_field_names = @_ ;
  
  $message = "undefined error" unless $message ;
  push @{$self->{_param_err_msgs}}, $message ;
  
  foreach (@param_field_names) {
    $self->{_param_err_names}{$_} = $message ;
  }

}

#######################################################
# Setup the job directory structure.
# We use tmpdir to create a unique directory name from
# the job name given by the user. This dir name then
# serves as the unique job_id.
# Argument:
# job_dir_server_path  abs path to the job storage dir
# Create:
# JOB_NAME_XXXXXXX/
#  data/
#  results/
#  status/
#  config/
#    tool_name (contains only name of this tool)
# 
#######################################################
sub setup_job_directory{
  my $self = shift ;
  my $job_storage_dir_server = shift || SrnaTools::Exception->throw('missing parameter job_dir_eerver') ;
  
  if (! -d $job_storage_dir_server || ! -w $job_storage_dir_server) {
    SrnaTools::Exception::FileAccess->throw(message=>'Job storage directoy on server is not a directory or not read/writable', log_msg=>"directoy: $job_storage_dir_server") ;
  }
  
  # Create the directory for the job
  # Name of directory follows the form
  # JobName_XXXXXXXXX (X=random)
  my $job_name = $self->job_name || 'MyJob' ;
  my $job_dir_template = $job_name.'_XXXXXXXXX' ;
  my $job_dir ;
  eval {
  $job_dir = tempdir( $job_dir_template, 
                      DIR => $job_storage_dir_server, 
                      CLEANUP => 0 ); # need to keep this
  # Make this world-writable
  chmod 0777 , $job_dir ;
  $self->job_working_dir($job_dir) ;
  } ;
  if ($@) {
    SrnaTools::Exception->throw(message => 'something went wrong trying to create a temporary directory in the job storage area', log_msg => $@) ;
  }
  # Now get the job_id, i.e. the name of
  # the newly created directoy
  $self->job_id( basename($job_dir) ) || SrnaTools::Exception->throw(message=>'Could not get job ID',log_msg=>"job_dir: $job_dir");
  
  # Create the subdirs for user data, results,
  # config files and status in the project dir
  foreach ( qw( data results status config lib) ) {
    my $dir = $job_dir.'/'.$_ ;
    mkdir $dir or SrnaTools::Exception::FileAccess->throw( "Could not create $dir sub-directory in job directory" );
    chmod 0777, $dir or SrnaTools::Exception::FileAccess->throw("Could not change permissions for $dir sub-directory in job directoy");
  }
  
  # Create a file with the name of the tool
  # in config dir
  open(NAME_FILE,'>',$job_dir.'/config/tool_name') or SrnaTools::Exception::FileAccess->throw( "Could not create a file in job config directory") ;
  print NAME_FILE $self->tool_name ;
  close NAME_FILE ;
  
  # create a shell script to submit remote
  # batch jobs
  if ($self->remote_mode) {
    my $script_file = $job_dir.'/lib/run_job.sh';
    open(SUBMISSION_SCRIPT, '>', $script_file);
    print SUBMISSION_SCRIPT $self->generate_submission_script ;
    close SUBMISSION_SCRIPT;
    chmod 0555, $script_file ;
  }

  return ;
} #setup_job_directory


#######################################################
# Delete the job working dir
#######################################################
sub delete{
  my $self = shift ;
  
  my $cmd = "rm -rf ".$self->job_working_dir;
  system($cmd);
  
} # delete


#######################################################
# Write user parameters into a file 'params' in the
# job's config directoy. Add job_name and id to params,
# this makes it easier to recreate the job from files
#######################################################
sub create_job_param_file{
  my $self = shift ;
  
  # we must have params at this stage
  if (! $self->params) {
    SrnaTools::Exception::JobPreparation->throw("Could not retrieve user-parameters for job when trying to write to params file") ;
  }
  
  my $config_dir = $self->job_working_dir.'/config' ;
  if ( ! -w  $config_dir ) {
    SrnaTools::Exception::FileAccess->throw(message=> "Can not write to config sub-directoy in job directory",log_msg=>"dir: $config_dir" );
  }
  my $params_file = $config_dir.'/params' ;
  
  $self->param('job_name',$self->job_name) ;
  $self->param('job_id',$self->job_id) ;
  
  $self->_write_hash_to_file($self->params, $params_file) ; 
  
} # create_job_param_file


#######################################################
# Create a copy of the tool config in a file in this
# job's working directoy. This is only used to process
# tools asynchronously so that they carry the config
# from the app that created them to wherever they are
# actually executed.
#######################################################
sub create_tool_config_file{
  my $self = shift ;
  my $tool_conf = $self->config ;
  
  if (! $self->config) {
    SrnaTools::Exception::JobPreparation->throw("Could not retrieve tool config for job when trying to write to config file") ;
  }
  my $config_dir = $self->job_working_dir.'/config' ;
  if ( ! -w  $config_dir ) {
    SrnaTools::Exception::FileAccess->throw(message=> "Can not write to config sub-directoy in job directory",log_msg=>"dir: $config_dir" );
  }
  my $config_file = $config_dir.'/tool.conf' ;
  
  $self->_write_hash_to_file($self->config, $config_file) ; 
  
} # create_tool_config_file


#######################################################
# Create a copy of the application config file for jobs
# that are processed remotely.
# Usage:
# $job->create_app_config_file( $config_hashref ) ;
#######################################################
sub create_app_config_file{
  my $self = shift ;
  my $config_hashref = shift ;
  
  if (! $config_hashref || ! ref $config_hashref eq 'HASH') {
    SrnaTools::Exception->throw(message => 'No hashref passed or not a hashref while creating copy of application conf file.', log_msg => "hashref: $config_hashref") ;
  }
  
  my $config_dir = $self->job_working_dir.'/config' ;
  if ( ! -w  $config_dir ) {
    SrnaTools::Exception::FileAccess->throw(message=> "Can not write to config sub-directoy in job directory",log_msg=>"dir: $config_dir" );
  }
  my $config_file = $config_dir.'/application.conf' ;
  
  $self->_write_hash_to_file($config_hashref, $config_file) ; 
  
} # create_app_config_file


#######################################################
# Use Data::Dumper to dump contents of a hashref into 
# the given file for configuration or user params. The
# format of the resulting file is Perl.
# NOTE: could use XML (e.g. XML::Simple) as alternative
# but not YAML because it is incompatible with the 
# Config::Auto parser 
# Usage:
# $job->write_to_config_file( $config_hashref, $file_path)
#######################################################
sub _write_hash_to_file{
  my ($self, $hashref, $file) = @_ ;
  if (! $hashref or ! ref $hashref eq 'HASH') {
    SrnaTools::Exception->throw(message=>"No hashref given or not a reference to a hash",log_msg=>"hashref: $hashref") ;
  }
  if (! $file) {
    SrnaTools::Exception->throw("No file path given") ;
  }
   
  my $output = Data::Dumper->Dump([$hashref]) ;
  open (FILE, '>', $file) or SrnaTools::Exception::FileAccess->throw(message=> "Problem writing config/parameter data to file: file is not writable", log_msg=>"file> $file");
  
  print FILE $output ;
  close FILE ;
  
} # write_to_config_file


#######################################################
# This method is to be implemented in those tools that
# need to display some result data in their results page
# beyond files in /results subdir. For example, the hp
# tool uses it to display a legend for the structure image
#######################################################
sub get_result_page_data{
  return undef;
}

#######################################################
# A shell script that can be submitted to the job queue
# to run the CliApp in batch mode on remote machine.
# This is also the place where parameters for the qsub
# command can be set
#######################################################
sub generate_submission_script{
  my $self = shift ;
  
  my $queue = $self->queue ;
  my $email = $self->app->admin_email;
  my $remote_root_dir = $self->app->root_dir($self->execution_env);
  my $remote_job_dir = $self->app->path_to('queue_job_dir',$self->execution_env);
  my $remote_path_job = $remote_job_dir.'/'.$self->job_id ;
  
  my $script =qq^#!/bin/csh
# Select queue
#\$ -q $queue.q

# Send mail if job aborted
# to webadmin address
#\$ -M $email
#\$ -m a

# NOTE this will have to be changed when 
# adapting to a different environment.
# All this is does is to make sure that
# all 3rd party binaries are in the PATH
module add srna-tools

# Run the job
$remote_root_dir/srna-tools.pl --run_batch_job $remote_path_job --queue_job_id \$JOB_ID --queue_hostname `hostname`
  ^;

  return $script ;
  
} # generate_submission_script


#######################################################
# A very basic email format validation
# Just checks for NAME@DOMAIN.SUFFIX
# We don't allow /:;* for security reasons
# even though they may be valid in emails
#######################################################
sub _email_is_valid{
  my ($self, $email) = @_ ;
  return 0 if $email=~m/[\/:;*]/ ;
  return ($email =~/^.+@.+\..+$/) ? 1 : 0 ;
}

#######################################################
# Convert this_string to ThisSting
# to guess Package/Module names
#######################################################
sub _camelcase{
  my ($self,$name) = @_ ;
  $name =~s/(\b|_)([a-z])/\u$2/ig ;
  return $name ;
}


# DO NOT DELETE
1 ;
