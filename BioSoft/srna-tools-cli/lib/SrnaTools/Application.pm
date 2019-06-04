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


######################################################################
#
# This is the Srna-Tools Application base class that handles common 
# functionality for the UEA sRNA toolkit, which is required by
# web and command line versions of the toolkit.
#
######################################################################

package SrnaTools::Application;
use strict;
use warnings;
use File::Temp qw( tempfile tempdir );
use YAML::Tiny ;
use Config::Auto ;
use Cwd 'abs_path';
use Exception::Class::TryCatch;


use SrnaTools::Exceptions;

# allowed environments (web, command line)
use constant ALLOWED_ENVIRONMENTS => qw( web cmd );

#######################################################
# Create a new app instance. Allowed params:
# mode : web/cmd
# cfg_file: path to app config file
#######################################################
sub new{
  my $class = shift ;
  my $self = {} ;
  bless $self, $class ;
  $self->_init(@_) ;
  return $self ;
} # new

#######################################################
# some additional setting up
#######################################################
sub _init{
  my $self = shift ;
  my %args = @_ ;
  $self->environment($args{environment}) if $args{environment} ;
  $self->cfg_file($args{cfg_file}) if $args{cfg_file};
  $self->{_tools_cfg_file} =  $args{tools_cfg_file} if $args{tools_cfg_file} ;
}

#######################################################
# Load config from a file. 
#######################################################
sub cfg_file{
  my ($self, $cfg_file) = @_ ;
  SrnaTools::Exception->throw("missing path to config file") unless $cfg_file ;
  SrnaTools::Exception->throw("config file $cfg_file is not readable") unless -r $cfg_file ;
  $self->{_cfg_file} = $cfg_file ;
  eval {
    $self->{_cfg} = Config::Auto::parse( $cfg_file, format => "perl" );
  } ;
  if ($@) {
  SrnaTools::Exception::Misconfiguration->throw(message=>"There has been a problem reading the config file",log_msg=> $@);
  }
}


######################################################
# Accessors
#######################################################
# Return the config hash
sub config{ $_[0]->{_cfg} || {} } 
sub cfg{ $_[0]->{_cfg} || {} } 
sub tools_cfg_file{ $_[0]->{_tools_cfg_file} } 

# An Application can handle a job
# and store the Job object
sub job{
  my ($self,$job) = @_ ;
  if ($job){
    if (! $job->isa('SrnaTools::Job') ) {
      SrnaTools::Exception->throw("Parameter is not a SrnaTools::Job object");
    }
    $self->{_job} = $job ;
  }
  $self->{_job};
}

#######################################################
# access application config data
#######################################################
sub tools{$_[0]->cfg->{tools} } # list fo available tools
sub debug{$_[0]->cfg->{debug} }
sub admin_email {$_[0]->cfg->{admin_email} }
sub max_number_jobs_queue {$_[0]->cfg->{max_number_jobs_queue} }
sub email_from {$_[0]->cfg->{email_from} }
sub no_job_submission {$_[0]->cfg->{no_job_submission} }
sub dont_delete_job_working_dir {$_[0]->cfg->{dont_delete_job_working_dir} }
sub root_dir{
  my ($self, $env) = @_ ;
  SrnaTools::Exception->throw("missing argument: environment (cluster/server)") unless $env ;
  return $self->cfg->{root_dir}{$env} ;
}
sub post_max{$_[0]->cfg->{post_max}}

sub queue_job_dir{
  my ($self, $env) = @_ ;
  SrnaTools::Exception->throw("missing argument: environment (cluster/server)") unless $env ;
  return $self->cfg->{queue_job_dir}{$env} ;
}

# Check with config that a given tool is available in
# this instance of the app.
# usage:
# $app->tool_is_available('mirprof') 
sub tool_is_available{
  my $self = shift ;
  my $tool_name = shift or SrnaTools::Exception->throw("missing tool name") ;
  my $cfg = $self->cfg ; # get config
  return grep ($_ eq $tool_name, @{$cfg->{tools}}) ? 1 : 0 ;
} # tool_is_available


# Retrieve path information from the config given a
# keyword for the path (e.g. 'queue_job_dir') and
# the name of the environment, e.g. cluster or server.
# register an exception if we can't retrieve it.
# Some paths need no environment, e.g. paths to web-only
# templates will only be on the server.
sub path_to{
  my $self = shift ;
  my $path_name = shift; # e.g. queue_job_dir
  my $env = shift || 'server' ; # cluster/server
  my $cfg = $self->cfg ; # get config
  
  SrnaTools::Exception->throw("No path name given") unless $path_name ;
  
  my $root_dir = $self->root_dir($env) ;
  SrnaTools::Exception->throw("No root directory given for $env in config") if !$root_dir;
  
  # some paths don't have cluster/server env
  # i.e. they are not references to hashes
  my $requested_path = ref  $cfg->{$path_name} eq 'HASH' ? $cfg->{$path_name}{$env} : $cfg->{$path_name} ;
  SrnaTools::Exception->throw(message=>"Problem with path configuration for $env",log_msg=>"Was asked for path '$path_name' in environment '$env', root_dir from config is '$root_dir'. Got this in config for this path: ". Data::Dumper->Dump([$cfg->{$path_name}])) if !$requested_path;
  return $root_dir.'/'.$requested_path ;

} # path to


# user_host_name is a string of the format
# user@host
# from the config file. We have one for server and
# one for the cluster
sub user_host_name{
  my $self = shift ;
  my $machine = shift ;
  SrnaTools::Exception->throw("no machine name given") unless $machine ;
  
  my $user_host_name = $self->cfg->{user_host_name}{$machine} ;
  if (! $user_host_name) {
     SrnaTools::Exception->throw("user_host_name for $machine not defined in config") ;
  } else {
    return $user_host_name;
  }
} # user_host_name


############################
# job dir usage allowance
# on server/cluster 
############################
sub usage_allowance{
  my $self = shift ;
  my $machine = shift ;
  SrnaTools::Exception->throw("machine name required") unless $machine ;
  return $self->cfg->{usage_allowance}{$machine} ;
} 


#######################################################
# read/set the environment (web or command line (cmd))
#######################################################
sub environment{
  my $self = shift ;
  my $arg = shift ;
  if ($arg) {
    if ( ! grep($_ eq $arg, ALLOWED_ENVIRONMENTS) ) {
      SrnaTools::Exception->throw($arg.' is not one of the allowed modes ('.join(',',ALLOWED_ENVIRONMENTS).').') ;
    } else {
      $self->{_environment} = $arg ;
    }
  } # if arg
  return  $self->{_environment} ;
} # mode


#######################################################
# Return a list of available tools as a hash ref with
# tool display names and short descriptions (if defined
# in the config file).
# This requires loading the Tool config for each tool in
# the app's tool_available list. For that we need the path
# to the tool config file. If this is not given, generate
# a list of tool names only.
# usage:
# foreach my $tool ( $app->tool_list('path/to/tools.config') ) {
#   print $tool->{display_name}.':  '. $tool->{description}."\n" ;
# }
######################################################
# TODO this should not do its own tool config parsing
#######################################################
sub tool_list{
  my $self = shift ;
  my $tools_config_file = shift || $self->tools_cfg_file ;
  my $cfg = $self->config ; # get config
  my $tool_cfg ;
  
  my @tools = () ;
  # Load the tool config file
  if (defined $tools_config_file) {
    SrnaTools::Exception::FileAccess->throw(message=> "could not find/read tool config file",log_msg=>"Trying to access '$tools_config_file'") unless -r $tools_config_file ;
    eval {
      $tool_cfg = Config::Auto::parse( $tools_config_file, format => "perl" );
    SrnaTools::Exception::Misconfiguration->throw(message=>"Problem reading tool config file",log_msg=>"$@") if $@ ; 
    } ;
    
  }
  
  # traverse the available tools and get
  # display name and description
  foreach my $tool ( @{$cfg->{tools}} ) {
    my $dname = $tool ;
    my $desc = 'no description available' ;
    if ($tool_cfg) {
      $dname = $tool_cfg->{$tool}{display_name} if $tool_cfg->{$tool}{display_name};
      $desc = $tool_cfg->{$tool}{description} if $tool_cfg->{$tool}{description};
    }
    push @tools, {cmd_name=> $tool, display_name => $dname, description => $desc} ;
  }

return @tools ;

} # tool_list


#######################################################
# As a generic way of handling SrnaTools::Exceptions,
# simply die of the error. Overwrite this method in 
# WebApp and CliApp
#######################################################
sub handle_error{
  my ($self, $err) = @_ ;
  if (!$err) {
    die "No error object passed" ;
  } else {
    $err->rethrow ; 
  }
}


#######################################################
# Parse error messages into an error name, which is
# used to display the right template, and some details
# (if given) that are adjusted according to the
# setting of 'debug' in the app config: in debug mode
# we add the log_msg to the details for output.
# We just fetch the actual class name of the error
# as its name, so:
# SrnaTools::Exception::Backend::RemoteHostNoResponse
# is turned into 'remote_host_no_response'
# Exceptions may have an additional field 'log_msg'
# that contains a string to be printed to the error
# log.
# Returns error_name, details, log_details
#######################################################
sub parse_error{
  my ($self, $err) = @_ ;
  my $err_name = 'undefined' ;
  my $details = '';
  my $log_details = '' ;
  
  if ($err->isa('SrnaTools::Exception') ){
    my $err_package = ref $err ;
    my ($class) = ($err_package=~/.*:{2}(\w+)$/) ;
    $err_name = $self->de_camelize($class) ;
    
    $details = $err->message if $err->message ;
    $log_details = $details ."; File: ".$err->file."; Line: ".$err->line."; Log-msg: ".($err->log_msg || '-');
    
    $details = $log_details if $self->debug;
   
   } else { # not an SrnaTools::Exception
      $details = $err if $self->debug ;
      $log_details = $err ;
   }
  
  return ($err_name, $details, $log_details) ;
} # parse_error


#######################################################
# Send a job to queue - this can be the local batch
# queue or the remote (cluster) job queue. This is
# decided by the tool config.
# tar|ssh (or cp for local batch queue)
# the whole directory to target. Not using
# any compression (tested with large sequence
# files and the compression actually slowed the
# process down) 
# Preserve permissions of files/dirs and
# compress with gzip when sending over network
# Use local cp if jobs are submitted to local queue.
# The tar command on the remote machine will return an error
# if the archive is incomplete. 
# Try 3 times before giving up.
# After a job is copied to the queue, we immediately
# any user data that was uploaded (to save space)
#######################################################
sub _send_job_to_queue_directory{
  my $self = shift ;
  my $job = $self->job ;
  
  SrnaTools::Exception->throw("Mo Job associated with App") unless $job;
  
  my $job_storage_dir = $self->path_to('job_storage_dir') ;
  my $job_dir_name = $job->job_id ;
  if (!$job_dir_name){
    SrnaTools::Exception->throw("Could not retrieve job ID and therefore could not find correct job directory");
  }
  my $env =  $job->execution_env ; # 'cluster' or 'server'
  my $target_dir = $self->path_to('queue_job_dir',$env) ;
  
  my $copy_cmd;
  my $set_remote_flag_cmd ; # for /status/awaiting_submission
  my $set_local_flag_cmd ; # for /status/sent_to_queue_directory
  
  $set_local_flag_cmd = 'touch '.$job_storage_dir.'/'.$job_dir_name.
                         '/status/sent_to_queue_directory';
  my $delete_upload_files_cmd = "rm $job_storage_dir/$job_dir_name/data/* 1>/dev/null 2>&1";
  
  # send tar archive to remote via ssh
  # $PIPESTATUS is non-zero if one of the
  # subcommands failed, otherwise this would
  # return zero if tar failed but ssh worked
  if ( $job->queue_remote_mode ) {
    $copy_cmd = "tar -p -c -C ".$job_storage_dir." -f - ".
            $job_dir_name ."| ssh ".$self->user_host_name($env).
            " tar xf - -C ".$target_dir . ' 1>/dev/null 2>&1 '.
            "; exit \$PIPESTATUS";
    $set_remote_flag_cmd = 'ssh '.$self->user_host_name($env ).
                    ' touch '.$target_dir.'/'.$job_dir_name.
                    '/status/awaiting_submission'.
                    ' 1>/dev/null 2>&1';
  } elsif ($job->queue_local_mode) {
    $copy_cmd = "cp -pr ". $job_storage_dir . '/' . $job_dir_name .
           " " . $target_dir ;
    $set_remote_flag_cmd = 'touch '.$target_dir.'/'.$job_dir_name.
                           '/status/awaiting_submission';
  } else {
    SrnaTools::Exception->throw("Was expecting a job in 'queue_local' or 'queue_remote' mode but job was neither of those.") ;
  }
  
  my $tried = 0 ;
  my $status = 1;
  my $errs ;
  while (1) {
    eval { $status = system($copy_cmd); };
    last if !$@ && $status == 0;
    if ( ++$tried == 3) {
      SrnaTools::Exception->throw(message=>"Failed to copy data to queue storage device: exceeded max. number of attempts. Please try again later.",log_msg=> "Tried 3 times to run this command: $copy_cmd") ;
    }
    sleep 5 ;
  }
  
  # Set the awaiting_submission flag
  system($set_remote_flag_cmd) == 0 or SrnaTools::Exception->throw(message=>"Failed to mark job as ready for submission to job queue. Please try again later.",log_msg=> "Failed to run this command: $set_remote_flag_cmd") ;
  
  system($set_local_flag_cmd) == 0 or SrnaTools::Exception->throw(message=>"Failed to mark job as 'sent to queue'. Please try again later.",log_msg=> "Failed to run this command: $set_local_flag_cmd") ;
  system($delete_upload_files_cmd); # not too critical if that failed now

} # _send_job_to_queue_directory


#######################################################
# Run a command on local or remote machine to submit
# a job to the queing system (batch or qsub)
#######################################################
sub _submit_job_to_queue{
  my $self = shift ;
  my $job = $self->job ;
  SrnaTools::Exception->throw("Mo Job associated with App") unless $job;
  
  if ($self->no_job_submission) {
    warn "!!! DEBUG option no_job_submission is active - job will be created and sent to remote host but not submitted to job queue !!!";
    return 0 ;
  }
  
  my $env =  $job->execution_env ; # 'cluster' or 'server'
  my $job_working_dir = $self->path_to('queue_job_dir',$env).'/'.$job->job_id ;
  my $cli_app = $self->root_dir($env).'/srna-tools.pl' ;
  
  my $submit_cmd;
  if ( $job->queue_remote_mode ) {
    $submit_cmd = "ssh ".$self->user_host_name($env)." 'source .bashrc; qsub $job_working_dir/lib/run_job.sh' 1>/dev/null 2>&1" ;
  } elsif ($job->queue_local_mode) {
    $submit_cmd = "echo \"$cli_app --run_batch_job $job_working_dir\" | batch 1>/dev/null 2>&1" ;
  } else {
    SrnaTools::Exception->throw("Was expecting a job in 'queue_local' or 'queue_remote' mode but job was neither of those.") ;
  }

  my $tried = 0 ;
  my $status = 1;
  my $errs ;
  while (1){
    eval { $status = system($submit_cmd); };
    last if !$@ && $status == 0;
    if ( ++$tried == 3) {
      SrnaTools::Exception->throw(message=>"Failed to submit job to queue: exceeded max. number of attempts. Please try again later.",log_msg=> "Tried 3 times to run this command: $submit_cmd, Error: $errs") ;
      # TODO delete the job remotely and write to error log
    }
    sleep 5 ;
  } 
  return (1-$status) ; # return true on success
  
} #_submit_job_to_queue

#######################################################
# Get a list of files in the job's result dir
#######################################################
sub get_job_result_files{
  my $self = shift ;
  my $job = $self->job ;
  SrnaTools::Exception->throw("Mo Job associated with App") unless $job;
  SrnaTools::Exception->throw("Could not retrieve job ID and therefore could not find correct job directory") unless $job->job_id;
  
  my $job_storage_dir = $self->path_to('job_storage_dir') ;
  my $job_results_dir = $job_storage_dir.'/'.$job->job_id .'/results';
  if (! -e $job_results_dir || ! -d $job_results_dir || ! -r $job_results_dir){
    SrnaTools::Exception->throw(message=>"The results directory for this job is not readable",log_msg=> "dir: $job_results_dir") ;
  }
  opendir(DIR, $job_results_dir);
  my @files = readdir(DIR);
  @files = grep ($_ !~/^\./, @files);
  closedir(DIR);
  
  return \@files ;
} # get_job_result_files

#######################################################
# Write job activity messages to log file
#######################################################
sub log_job_activity{
   my $self = shift ;
   my $msg = shift || 'undefined activity - got no message' ;
   my $job = $self->job;
   SrnaTools::Exception->throw("App is not associated with a job") unless $job ;
   
   my $queue_job_id = $self->param('queue_job_id') || 'n/a';
   my $queue_hostname = $self->param('queue_hostname') || 'localhost';
   my $log_file = $self->path_to('log_file', $job->execution_env) ;
   my $time = `date` || 'could not get time' ;
   my $tool_name = $job->tool_name ;
   my $job_id = $job->job_id ;
   
   my $cmd =qq^echo "`date`: tool: $tool_name, job_id: $job_id, queue-job-ID: $queue_job_id, host: $queue_hostname -> $msg" >> $log_file^;

   system($cmd) ;
   
} # log_job_activity


#######################################################
# Get current status of job from status file. This is
# also called by a AJAX snippet in the page for updates
# This returns 'completed' if the job is completed,
# regardless of errors. Use _get_last_job_status
# for fetching the last recorded status from
# status_messages (never 'completed')
#######################################################
sub _get_job_status{
  my $self = shift ;
  my $job_id = shift || $self->query->param('jobid');
  my $status ;

  if (! $job_id) {
    SrnaTools::Exception->throw("subroutine called with missing parameter job_id");
  }
  
  # shortcut - if this file exists
  # the job is completed
  my $completed_file = $self->path_to('job_storage_dir').'/'.$job_id.'/status/completed' ;
  if (-e $completed_file) {
    return "completed" ;
  }

  $status = $self->_get_last_status_msg ;
  $status = 'queued...' if $status eq 'none' or !$status ;
  
  return $status  ;
} # _get_job_status


######################################################
# Get current status of job from status file. This is
# also called by a AJAX snippet in the page for updates
# This returns 'completed' if the job is completed,
# regardless of errors. Use _get_last_status_msg
# for fetching the last recorded status from
# status_messages (never 'completed')
#######################################################
sub _get_last_status_msg{
  my $self = shift ;
  my $job_id = shift || $self->query->param('jobid');
  my $status = 'none' ;

  if (! $job_id) {
    SrnaTools::Exception->throw("subroutine called with missing parameter job_id");
  }
  
  my $status_file = $self->path_to('job_storage_dir').'/'.$job_id.'/status/status_messages' ;
  ($status_file)=($status_file=~/(.+)/) ;
  
  # Get the last line from this file
  if (-r $status_file) {
    my $cmd = "tail -n 1 $status_file" ;
    eval {$status = `$cmd`;};
    if ($@){
      SrnaTools::Exception->throw(message=>"something went wrong trying to read the status of this job",log_msg=>"The following command failed although the file was readable: $cmd. Error: $@")
    }
  }
  return $status ;
} # _get_job_status


#######################################################
# Transfer a batch job back to the 'server'.
# The transfer protocol depends on the mode of the job,
# with local jobs using simple copy and remote jobs
# tar and ssh. 'server' means the machine on which the
# job was generated - for local instant or batch queue
# jobs this is the same machine as where they are
# executed. This method is CliApp specific because batch
# jobs are only handled by the CliApp.
#######################################################
sub _transfer_job_to_server{
  my $self = shift ;
  my $job = $self->job ;
  if (!$job ) {
    SrnaTools::Exception->throw("App is not associated with a job object") ;
  }
  
  my $working_dir = $job->job_working_dir ;
  my $job_storage_dir_server = $self->path_to('job_storage_dir');
  my $target_dir = $job_storage_dir_server.'/'.$job->job_id ;
  
  my $copy_cmd;
  my $set_remote_flag_cmd ; # for /status/completed
  my $delete_cmd = "rm -rf $working_dir";
  
  if ( $job->queue_remote_mode ) {
    $copy_cmd = "tar -p -c -C ".$working_dir." -f - results ".
            "| ssh ".$self->user_host_name('server').
            " tar xf - -C ".$target_dir . ' 1>/dev/null 2>&1 '.
            "; exit \$PIPESTATUS";
    $set_remote_flag_cmd = 'ssh '.$self->user_host_name('server').
                    ' touch '.$target_dir.'/status/completed'.
                    ' 1>/dev/null 2>&1';
  } elsif ($job->queue_local_mode) {
    $copy_cmd = "cp -r ". $working_dir . '/results ' . $target_dir. ' 1>/dev/null 2>&1 ' ;
    $set_remote_flag_cmd = 'touch '.$target_dir.'/status/completed';
  } else {
    SrnaTools::Exception::ModuleExecution->throw("Was expecting a job in 'queue_local' or 'queue_remote' mode but job was neither of those.") ;
  }
  
  my $tried = 0 ;
  my $status = 1;
  my $errs ;
  do {
    eval { $status = system($copy_cmd); };
    $errs = $@ ;
    if (++$tried == 3) {
      # can't copy to server, so no way to inform
      # running process. Just try to log this
      # locally and die
      my $msg = "Failed to copy data to server: exceeded max. number of attempts. Command: $copy_cmd" ;
      eval{ $self->log_job_activity("ERROR: $msg") }; 
      die $msg ;
    }
    sleep 5 ;
  } until ($status == 0 && !$errs );
  
  # Set the awaiting_submission flag
  system($set_remote_flag_cmd) == 0 or SrnaTools::Exception->throw(message=>"Failed to mark job as completed on server.",log_msg=> "Failed to run this command: $set_remote_flag_cmd") ;
  unless ($self->dont_delete_job_working_dir){
    system($delete_cmd)==0 or warn "Could not delete completd job $working_dir" ;
  }
  
}


#######################################################
# Update status of a running job in the storage directory
# on the server. This can be called from back-end modules
# Usage:
# $app->update_job_status(Job_object, message).
# Failure shouldn't throw an error because it is not
# ciritcal.
#######################################################
sub update_job_status{
  my $self = shift ;
  my $job = shift ;
  my $status_msg = shift || '';
  my $target_file = $job->job_id.'/status/status_messages' ;
  
  $self->send_msg_to_file_on_server($job, $status_msg, $target_file) ;
}

#######################################################
# Append a message to a file in the job dir on the
# machine where the job originated. This is used to
# update the status and log errors of jobs that run
# in batch mode on remote or local host.
# Strip message of all quotes - they can cause problems
# because we need to append the msg with 'echo' to a 
# file and nested quotes can cause this to fail
#######################################################
sub send_msg_to_file_on_server{
  my $self = shift ;
  my $job = shift ;
  my $msg = shift ;
  my $target_file = shift ;
  
  # sanitize - replace quotes with _
  # and make sure we have no newlines
  # and trailing whitespace
  $msg=~s/["'`]/_/g ;
  $msg=~s/\n/ /g ;
  $msg=~s/\s+$// ;
  
  if (!$job || ! $job->isa('SrnaTools::Job') ) {
    SrnaTools::Exception->throw("Parameter 'job' missing or not a SrnaTools::Job object") ;
  }
  
  if (!$msg || !$target_file) {
    SrnaTools::Exception->throw("Missing parameters: need 'msg' and 'target_file'");
  }
  
  my $job_storage_dir_server = $self->path_to('job_storage_dir');
  $target_file = $job_storage_dir_server .'/'. $target_file ;
  my $cmd ;
  
  if ( $job->queue_remote_mode ) {
    $cmd = "ssh ".$self->user_host_name('server').
           " 'echo \"$msg\" >> $target_file' 1>/dev/null 2>&1";
  } elsif ($job->queue_local_mode) {
    $cmd = "echo '$msg' >> $target_file 2>/dev/null";
  } elsif ($job->instant_mode && $self->environment eq 'cmd') {
    warn "$msg\n";
    return ;
  }
  
  my $tried = 0 ;
  my $status = 1;
  my $errs ;
  do {
    eval { $status = system($cmd); };
    $errs = $@ ;
    if (++$tried == 5) {
      # Pointless to raise another proper error because we can't communicate
      # with the server anymore
      warn "Seem to have lost connection between machines - can't update job status from remote batch (tried 5 times). Command:".$cmd;
      exit 1;
    }
    sleep 3 * $tried + 1 ;
  } until ( $status == 0 && !$errs);
  
  return (1 - $status ); # true on success of cmd
}


#######################################################
# Batch jobs log their runtime errors to a file
# in the job directory on the host where the job
# was generated. These can then be shown e.g. in the
# status page of the web site.
# The log file is JOBDIR/status/errors
#######################################################
sub handle_queue_job_error{
  my $self = shift ;
  my $err = shift || 'unknown exception';
  my $job = $self->job ;
  
  # TODO an error here would just result in a 'die'
  # how to handle this? Maybe try to log it at least
  SrnaTools::Exception->throw("App is not associated with a job") unless $job;
  my $target_file = $job->job_id.'/status/errors' ;
  $self->send_msg_to_file_on_server($job, $err, $target_file) ;
  $self->_transfer_job_to_server ; 
  if ($job->param('email') ) {
    $self->_send_completion_notification_from_server;
  }
}


#######################################################
# Display all errors in a similar sub-layout into which
# we render the error page specified here.
# The name of the error is passed as param 'error' and
# there must be a template file. Naming conventions
# for the error html partial:
# NameOfError_error.tt
# for example to set up an out-of memory error:
# URL param:
# error=out_of_memory
# create a template in html_templates and cli_templates
# out_of_memory_error.tt.
#
# This method is used by Web and Cli implementations of
# the App - the only difference is that the CliApp uses
# templates from the cli_templates dir instead of
# html_templates, which is configured on setup of the App
#######################################################
sub error_page{
  my $self = shift;
  my $error_name = shift;
  my $details = shift ;
  my $params = {} ;
  
  # Select an error html partial or use a generic one
  $params->{partial} = $error_name ? $error_name.'_error.tt' : 'undefined_error.tt' ;  
  
  # Reduce verbosity of "details" unless
  # in debug mode
  if ( $details ) {
    chomp $details ;
    if (! $self->debug) { 
      $details=~s/ at.+line \d+// ;
      $details=~s/\.*\s*$/./ ; # just one full stop and no trailing spaces
    }
    $params->{details} = $details ;
  }
  
  # Use the error_page layout for this. 
  # Retry with generic template if the 
  # specific template wasn't found
  my $output ;
  eval {
    $output = $self->tt_process('error_page.tt', $params );
  } ;
  if ($@) {
   if ($@=~/file error .+ not found/) {
    my $unrecogn_type = $params->{error} ;
    $params->{partial} = 'undefined_error.tt' ;
    $params->{details} .= " The error type that was returned ('$unrecogn_type') is not recognised";
    $output = $self->tt_process('error_page.tt', $params ); # try again
    } else { # can't even find 'undefined_error.tt' template
      die "Can't find error page template: ". $@;
    }
  }
  return $output ;
} # error_page


#######################################################
# Call the CliApp on the server and let it send an
# email to the user. This is not critical so we just
# return if it didn't work
#######################################################
sub _send_completion_notification_from_server{
  my $self = shift ;
  my $job = $self->job ;
  
  eval {
  return if !$job;
  
  my $job_storage_dir_server = $self->path_to('job_storage_dir');
  my $target_dir = $job_storage_dir_server.'/'.$job->job_id ;
  my $cli_app = $self->root_dir('server').'/srna-tools.pl' ;
  
  my $cmd ;
  
  if ( $job->queue_remote_mode ) {
    $cmd = "ssh ".$self->user_host_name('server').
            " $cli_app --send_email_notification $target_dir" ;
  } else {
    $cmd = " $cli_app --send_email_notification $target_dir" ;
  } 
  #$self->log_job_activity("send email: $cmd");
  my $tried = 0 ;
  my $status = 1;
  do {
    eval { $status = system($cmd); };
    last if ++$tried == 3;
    sleep 5 ;
  } until ($status == 0 );
  
  }; # eval
}


#######################################################
# Send email notification for completed jobs
#######################################################
sub do_send_completion_notification{
  my $self = shift;
  my $job_dir = shift;
  
  return unless $job_dir;
  
  # this action is not critical and should not lead
  # to exceptions if it fails
  eval{
  use Mail::Sendmail ;
  
  my $job = SrnaTools::JobFactory->create_job_from_config(
    app             => $self,
    job_dir         => $job_dir,
    skip_app_conf   => 1, 
  ) ;
  return unless $job;
  
  my $to = $job->param('email');
  return unless $to;
  
  my $tool_name = $job->display_name;
  my $job_name = $job->job_name;
  my $link = $job->param('result_page_link') || '## ERROR generating link ##' ;
  my $from = $self->email_from ;
  
  my $subject = "UEA sRNA-tools server: job $job_name completed" ;
  my $body = <<BODY ;
Dear user of the UEA sRNA tools,

The following job has been completed:
Job  : $job_name
Tool : $tool_name

Please follow the following link to retrieve the results:
$link
        
Please ensure that the entire link appears in the address bar of your browser. If that is not the case, please copy and paste the full link into the address bar manually.
  
BODY

  my %mail = (
    From    => $from,
    To      => $to,
    Subject => $subject,
    Message => $body,
  );
  sendmail(%mail);
  
  }; # eval
} # send_completion_notification


#######################################################
# Send an email to the admins
#######################################################
sub send_email_to_admin{
  my ($self, $msg) = @_;
  $msg = 'instruction received to send an email to admins but no message was given - D\'OH' unless $msg;
  
  eval{
  use Mail::Sendmail ;

  my $to = $self->admin_email;
  return unless $to;
  
  my $from = $self->email_from ;
  
  my $subject = "UEA sRNA-tools server: exception/warning" ;
  my $body = <<BODY ;

The following warning/exception message was received:

$msg

BODY

  my %mail = (
    From    => $from,
    To      => $to,
    Subject => $subject,
    Message => $body,
  );
  sendmail(%mail);
  
  }; # eval
  
}

#######################################################
# Turn 'MyClassName' into 'my_class_name'
#######################################################
sub de_camelize{
  my ($self, $string) = @_ ;
  return unless $string ;
  $string=~s/([A-Z])/_\l$1/g ;
  $string=~s/^_// ;
  return $string ;
}


# DO NOT DELETE
1;