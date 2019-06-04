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


####################################################
#
# This package contains the web-app for the
# srna-tools websites.
# It uses CGI::Application as a web framework and
# also inherits from SrnaTools::Application
#
####################################################

package SrnaTools::Application::WebApp;
use strict;
use base qw( CGI::Application SrnaTools::Application ); # don't change order!
use CGI::Application::Plugin::TT;
use CGI::Application::Plugin::Forward;
use CGI::Application::Plugin::Redirect;
use CGI::Application::Plugin::FillInForm (qw/fill_form/);

use Exception::Class::TryCatch;
use File::Copy ; 
use File::stat ; # only for getting job submission time
use Cwd 'abs_path';

use lib "../lib";
use SrnaTools::JobFactory ;

# config files
use constant APP_CONFIG_FILE   => '../config/application.conf' ;
use constant TOOLS_CONFIG_FILE => '../config/tools.conf' ;

# For links: this is the relative path to the
# cgi script that drives the site from itself
use constant THIS_SCRIPT_LINK => '../cgi-bin/srna-tools.cgi' ;
use constant THIS_SCRIPT_ABS_LINK => 'http://'.$ENV{SERVER_NAME}.'/'.$ENV{SCRIPT_NAME};

use Data::Dumper ;
## Use this for debugging:
# Put this anywhere into the main eval block
# of a runmode to see the data stored
# in an object (e.g. $job):
# $html= "<pre>".Data::Dumper->Dump([$job])."</pre>" ;return; 


######################################################
# Define our run modes
######################################################
sub setup {
    my $self = shift;
    $self->start_mode('document');
    $self->error_mode('error');
    $self->mode_param('rm');
    $self->run_modes(
      'input_form'   => 'input_form',
      'prepare_job'  => 'prepare_job',
      'submit_job'   => 'submit_job',
      'results'      => 'results',     # display job results
      'document'     => 'static_page', # render document named in parameter
      'error_page'   => 'error_page',  # render error_page named in parameter
      'job_status'   => '_get_job_status', # for AJAX calls
    );
}

#######################################################
# Execute the following before we execute the 
# requested run mode
#######################################################
sub cgiapp_prerun {
  my $self = shift;
  
  eval {
  # Read the config file
  # TODO handle config parse error better.
  # Currently, this is only logged to error-log
  # but the message on screen is about the missing
  # template file
  $self->cfg_file(APP_CONFIG_FILE) ;
  $self->environment('web') ;
  
  # untaint environment variables 
  $ENV{PATH} = "/bin:/usr/bin:/usr/local/bin" ;
  
  # configure Template::Toolkit
  $self->tt_config(
    TEMPLATE_OPTIONS => {
      INCLUDE_PATH => [
        $self->path_to('html_template_path'),
        $self->path_to('html_template_path').'/core',
        $self->path_to('html_template_path').'/site_spec_content',
        $self->path_to('html_template_path').'/core/errors',
        $self->path_to('html_template_path').'/core/layouts',
      ],
      ABSOLUTE     => 1,
      DEBUG => 1
    }
  ) ;
  
  # Add some standard params to all templates
  $self->tt_params(
    title           => $self->title,
    sidenav_partial => $self->sidenav_partial,
    static_doc_link => \&static_doc_link,
    image_src       => \&image_src,
    question_mark_link => \&question_mark_link,
    files_link      => \&files_link,
    tool_link       => \&tool_link,
  );

  } ;
  # catch errors
  if ( catch my $err ) {
  #return "<pre>".Data::Dumper->Dump([$err])."</pre>" ;
    $self->handle_error($err);
  }
} #cgiapp_prerun



#######################################################
# Process any fatal errors.
# NOTE this is only called if the script dies, e.g.
# as a result of a syntax error. For all other cases
# there are proper error pages that are rendered inside
# the general layout using the 'error_page' run mode.
# Only drop a short message into the browser but log
# the full error message.
#######################################################
sub error {
    my $self  = shift;
    my $error = shift;
    
    # don't be verbose in web page but
    # write the full message into the logs
    my $display_msg = ($error=~/syntax/)   ? 'Syntax error' :
                      ($error=~/file error/) ? 'file error: could be a missing template' : 
                      'Run time error' ;
    warn $error ;
    
    return "<html><body><div id=\"errorBox\" ><h2>There has been an error</h2><p> $display_msg. Details have been written to error log</p><p>please contact webmaster</p></div></body></html>";
} #error


#######################################################
# Display static content. By default, this will be the
# the home page. The page name to display is passed in
# as the page argument (without the suffix), so for the
# about page, the cgi parameter would be "page=about"
# and we would then need a file called 'about.tt' in
# our templates folder.
# TODO: the static templates should have their own folder
# which should be given in the config (or use some default)
# and added to the file path here
#######################################################
sub static_page{
  my $self = shift;
  my $errs = shift;
  my $q    = $self->query;
  
  my $params = {
    page => $q->param('page') || 'start_page',
  } ; # put additional page params here
  
  $params->{page}.='.tt'; 
  
  my $html = $self->tt_process('static_page.tt', $params ) || $self->handle_error(SrnaTools::Exception::TemplateError->throw("Problem processing static_page.tt".$self->tt_obj->error() ));
  return $html;
  
} # static_page


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
# create a template in html_templates
# out_of_memory_error.tt.
#
# This is handled by the Application class.
#######################################################
sub error_page{
  my $self = shift ;
  $self->SUPER::error_page(@_);
} # error_page


#######################################################
# Display job input forms.
# Before we can accept jobs we need to make a few
# health checks on the server and cluster to ensure
# that both are up and storage not overflowing. Also
# check for scheduled down time and stop accepting
# jobs if we are within a down time period.
#######################################################
sub input_form{
  my $self = shift;
  my $err_msgs = shift ;
  my $err_names = shift ;
  my $q    = $self->query;
  my $final_html ;
  
  eval {
  my $tool_name = $q->param('tool') ;
  if (! $self->tool_is_available($tool_name)) {
    SrnaTools::Exception::ToolNotAvailable->throw;
  }
  
  # Initialise a Job object for the requested tool
  # to access tool config for input field options
  my $job = SrnaTools::JobFactory->create_job(
    app            => $self,
    tool_name      => $tool_name, 
    tool_conf_file => TOOLS_CONFIG_FILE, 
  ) ;
  $self->job($job) ; # store connection to job
  
  $self->_run_site_health_checks($job->execution_env) unless $job->instant_mode;
  
  my $form_partial = $job->form_partial ;
  my $head_partial = $job->head_partial ;
  
  if (!$form_partial || !$head_partial) {
    SrnaTools::Exception::Misconfiguration->throw("Form partial $form_partial or header partial $head_partial not found as specified in config") ;
  }
  
  # Some parameters for the html template
  my $page_params = {
    head_partial    => $head_partial,
    content_partial => $form_partial,
    error_messages  => $err_msgs, 
    error_names     => $err_names, # names of fields that produced errors
    tool_config     => $job->cfg, 
    app_config      => $self->cfg,
  } ;
   
  # Generate the page and fill in forms
  my $html = $self->tt_process('tool_page.tt', $page_params ) || SrnaTools::Exception::TemplateError->throw("Problem processing tool_page.tt".$self->tt_obj->error() );
  $final_html = $self->fill_form( \$html, $q ) || SrnaTools::Exception::TemplateError->throw( "Failed to auto-fill input form");
  } ; #eval
  
  # catch errors or render page
  if ( catch my $err ) {
    $self->handle_error($err);
  } else {
    return $final_html ;
  }
} # input_form


#######################################################
# This run mode is called from the input forms. It 
# takes care of evaluating the provided parameters and
# either returns the user to the input form with errors
# if there are any, presents the actual submission page
# or a page that asks the user to wait for a link sent
# by email (for submitting the job)
#######################################################
sub prepare_job{
  my $self = shift;
  my $errs = shift;
  my $q    = $self->query;
  my $html ;
  my $job ;
  eval {
  my $tool_name = $q->param('tool') ;
  if (! $self->tool_is_available($tool_name)) {
    SrnaTools::Exception::ToolNotAvailable->throw;
  }
  
  my $user_params = $self->_parse_url_params;
  
  # Set some aditional params
  # that are always needed in WebApp
  $user_params->{do_uncompression} = 1 ;
  
  $job = SrnaTools::JobFactory->create_job(
    app            => $self,
    tool_name      => $tool_name, 
    tool_conf_file => TOOLS_CONFIG_FILE,
    environment    => $self->environment,
    params         => $user_params,
  ) ;
  $self->job($job) ; # store connection to job
  
  $self->_run_site_health_checks($job->execution_env) unless $job->instant_mode;
  
  my $post_max_mb = $self->post_max / 1e6 ;
  my $uploaded_mb = ($ENV{'CONTENT_LENGTH'} / 1048576 ) ;
  if ($uploaded_mb > $post_max_mb){
    my $post_max_str = sprintf("%.1f", $post_max_mb) ;
    my $uploaded_str = sprintf("%.4f", $uploaded_mb );
    SrnaTools::Exception::Fileupload->throw("trying to upload $uploaded_str MB. The maximum upload size is $post_max_str MB");
  }
  
  # validate the submitted parameters.
  # re-display input form with errors if neccessary
  if ( ! $job->validate_user_parameters) {
    my $err_msgs = $job->param_err_msgs; # error messages
    my $err_names = $job->param_err_names; # input field names for highlighting
    $html = $self->forward('input_form', $err_msgs, $err_names) ;
    return;
  }
  
  $job->setup_job_directory( $self->path_to('job_storage_dir') );
  
  my $results_link = THIS_SCRIPT_ABS_LINK.'?rm=results&jobid='.$job->job_id ;
  $job->param('result_page_link', $results_link);
  
  # Let job create local copies of
  # config files if not in 'instant' mode
#  if ( ! $job->instant_mode ) { 
    $job->create_job_param_file ;
    $job->create_tool_config_file ;
    $job->create_app_config_file($self->config) ;
#  }
  
  # Handle data files if any (configured for each
  # tool in config items data_file_params)
  # This is a CGI based app, so all files must
  # be uploaded to local folder. 
  if (ref $job->data_file_params eq 'ARRAY') {
    foreach my $file ( @{$job->data_file_params} ) {
      $self->_upload_file($file, $job->job_working_dir.'/data')
    }
  }
  
   if ($job->requires_email) {
      my $page_params = {
        head_partial    => $job->head_partial,
        content_partial => 'job_prepared_await_email.tt',
        job_name        => $job->job_name,
        email_addr      => $job->param('email'),
      } ;
      $html = $self->tt_process('tool_page.tt', $page_params ) || SrnaTools::Exception::TemplateError->throw("Problem processing tool_page.tt".$self->tt_obj->error() );
    } elsif ($job->instant_mode){ # go straight to run_job
      $job->run ;
      $self->redirect(THIS_SCRIPT_LINK.'?rm=results&jobid='.$job->job_id) ;
    } else { # submit and track the job
      $self->redirect(THIS_SCRIPT_LINK.'?rm=submit_job&jobid='.$job->job_id) ;
    }
    
  } ; # eval
   
  # catch errors or render page
  if ( catch my $err ) {
    $self->handle_error($err);
  } else {
    return $html ;
  }
} # prepare_job


#######################################################
# This run method submits the jobs and also updates
# its status until results are ready. Only asynchronous
# jobs need to use this.
#######################################################
sub submit_job{
  my $self = shift;
  my $job_id = $self->query->param('jobid') ;
  my $html ;
  
  eval {
  # must have a job id
  if (! $job_id ) {
    SrnaTools::Exception::Url::UrlParameterMissing->throw(message=>"job ID is required",log_msg=>"parameter 'jobid' missing or empty - jobid: $job_id") ;
  }
  
  $self->_job_dir_health_check($job_id);
  
  my $job = SrnaTools::JobFactory->create_job_from_config(
    app             => $self,
    job_dir         => $self->path_to('job_storage_dir').'/'.$job_id,
    skip_app_conf   => 1, 
  ) ;
  $self->job($job) ; # store connection to job
  
  # links to results page
  # and job_status (for AJAX calls)
  my $status_update_link = THIS_SCRIPT_LINK.'?rm=job_status&jobid='.$job_id ;
  my $results_link = THIS_SCRIPT_LINK.'?rm=results&jobid='.$job_id ;
  $job->param('result_page_link', $results_link);
  
  # Have we sent this to the queue? If not do it now
  # Also send the command to enter the job into the queue
  if (! $self->_job_sent_to_queue_directory($job_id)) {
    $self->_send_job_to_queue_directory ;
    $self->_submit_job_to_queue ;
  }

  # Get the header for this tool (from config)  
  my $head_partial = $job->head_partial ;
  if (!$head_partial) {
    SrnaTools::Exception::Misconfiguration->throw("header partial $head_partial not found as specified in config") ;
  }
  
  my $job_status = $self->_get_job_status($job_id) ;
  my $time_started = $self->_get_job_submission_time($job_id) ;
  
  my $result_link_visible = $job_status eq 'completed' ? 1 : 0 ;
  my $page_params = {
    head_partial          => $head_partial,
    content_partial       => 'job_status_page.tt',
    result_link_visible   => $result_link_visible,
    job_name              => $job->job_name,
    job_id                => $job_id,
    job_status            => $job_status,
    job_submission_time   => $time_started,
    status_update_link    => $status_update_link,
    results_link          => $results_link,
  } ;
  $html = $self->tt_process('tool_page.tt', $page_params ) || SrnaTools::Exception::TemplateError->throw("Problem processing tool_page.tt".$self->tt_obj->error() );

  } ; # eval
  
  # catch errors or render page
  if ( catch my $err ) {
    $self->handle_error($err);
  } else {
    return $html ;
  } 
  
} #run_job


#######################################################
# This run method displays job results. 
#######################################################
sub results{
  my $self = shift;
  my $job_id = $self->query->param('jobid') ;
  my $html ;
  
  eval {
  # must have a job id
  if (! $job_id ) {
    SrnaTools::Exception::Url::UrlParameterMissing->throw(message=>"job ID is required",log_msg=>"parameter 'jobid' missing or empty - jobid: $job_id") ;
  }
  
  # run some error checks on the job directory
  $self->_job_dir_health_check($job_id);
  
  # Recreate the job object from config data
  my $job = SrnaTools::JobFactory->create_job_from_config(
    app             => $self,
    job_dir         => $self->path_to('job_storage_dir').'/'.$job_id,
    skip_app_conf   => 1, 
  ) ;
  $self->job($job) ; # store connection to job
  
  # job must be completed, otherwise we need
  # to go back to the submit page
  # Instant jobs are always considered complete
  if (!$job->instant_mode && !$self->_job_complete($job_id) ) {
    $self->redirect(THIS_SCRIPT_LINK.'?rm=submit_job&jobid='.$job_id) ;
  }
  
  # Get the header for this tool (from config)  
  my $head_partial = $job->head_partial ;
  if (!$head_partial) {
    SrnaTools::Exception::Misconfiguration->throw("header partial $head_partial not found as specified in config") ;
  }
  
  my $err_msg = $self->_job_error_msg ;
  
  if ($err_msg) {
    my $last_status =  $self->_get_last_status_msg($job_id);
    my $page_params = {
      head_partial  => $head_partial,
      content_partial => 'job_error_page.tt',
      error_msg   => $err_msg,
      last_status => $last_status,
    } ;
    $html = $self->tt_process('tool_page.tt', $page_params ) || SrnaTools::Exception::TemplateError->throw("Problem processing tool_page.tt".$self->tt_obj->error() );
  
  } else { # no errors generate results page
    my $result_files = $self->get_job_result_files ;
    my $result_page_partial = $job->result_page_partial || 'generic_results_page.tt' ;
    my $job_results_link = '../jobs/'.$job_id.'/results/' ;
    my $tool_name = $job->tool_name;
    my $result_page_data = $job->get_result_page_data ;
    
    
    my $page_params = {
      head_partial     => $head_partial,
      tool_name        => $tool_name,
      content_partial  => $result_page_partial,
      result_files     => $result_files,
      job_results_link => $job_results_link,
      result_page_data => $result_page_data,
    } ;
    
    $html = $self->tt_process('tool_page.tt', $page_params ) || SrnaTools::Exception::TemplateError->throw("Problem processing tool_page.tt".$self->tt_obj->error() );
  }
  
  } ; # eval
  
  # catch errors or render page
  if ( catch my $err ) {
    $self->handle_error($err);
  } else {
    return $html ;
  } 
  
} #results


#######################################################
# Run some health checks on the hardware and the site
# in general. This includes checking disk space quotas
# and downtime schedules. 
# Return hash of error name and (optional) details
# that can be used to forward the run mode to 'error_page'
# Mode is 
#######################################################
sub _run_site_health_checks{
  my $self = shift ;
  my $execution_env = shift ; # cluster, server, local_queue
  SrnaTools::Exception->throw('no mode parameter passed to subroutine') if !$execution_env ;
  
  # check if we are in a scheduled downtime
  if (my $downtime_end = $self->_within_scheduled_downtime ){
    SrnaTools::Exception::Downtime->throw($downtime_end);
  }
  
  # Check disk quota and response from cluster backend
  my $disk_usage = $self->_get_job_dir_usage($execution_env) ;
  if (!defined $disk_usage) { # remote host didn't respond
    SrnaTools::Exception::Backend::RemoteHostNoResponse->throw;
  } else {
    my $allowance = $self->usage_allowance($execution_env) ;
    if (! defined $allowance) {
      SrnaTools::Exception::Misconfiguration->throw("usage allowance not defined for excetuion environment '$execution_env' in config");
    }
    
    if ($disk_usage > $allowance) {
      SrnaTools::Exception::ExceedJobDirQuota->throw;
    }
  }
  
  my $max_number_jobs_in_queue = $self->max_number_jobs_queue || 100;
  if ($execution_env eq 'cluster' && $self->_number_jobs_in_remote_queue >= $max_number_jobs_in_queue ){
    $self->send_email_to_admin("A job was rejected because the maximum number of jobs in queue ($max_number_jobs_in_queue) is exceeded on cluster. Please check the job queues on the cluster.");
    SrnaTools::Exception::ExceedJobNumberQueue->throw;
  }
  return 1 ; # all clear
} # _run_site_health_checks



#######################################################
# Run some checks on the job directory to make sure it
# is complete. We need a job id, a directory with that name
# must exist in the server job-directory and it must be 
# readable. The job directory itself must be available
# and we must have a results dir and a config dir
# inside that must be readable.
#######################################################
sub _job_dir_health_check{
  my $self = shift ;
  my $job_id = shift ;
  my $job_storage_dir = $self->path_to('job_storage_dir') ;
  
  if (! $job_id) {
    SrnaTools::Exception->throw("subroutine called with missing parameter job_id");
  }
 
  # prevent users from accessing random files
  my $sanitized_job_id = $self->_sanitize_param($job_id) ;
  if ( not $job_id eq $sanitized_job_id ) {
    SrnaTools::Exception::Url->throw(message=>"job_id parameter is incomplete or may have been modified.",log_msg=>"parameter 'jobid' did not pass the security check. This may indicate that a user has modified the parameter value that was given in a link and it may include unsafe characters that would allow access to files outside the job dir. The parameter value for job_id was: $job_id, sanitized_job_id: $sanitized_job_id") ;
  } 
  
  
  if (! $job_storage_dir) {
    SrnaTools::Exception::Misconfiguration->throw("path to job storage directory could not be found in application config");
  }
  
  if (! -d $job_storage_dir || ! -r $job_storage_dir) {
    SrnaTools::Exception::FileAccess->throw(-message=>"job storage directory on server is not accessible",log_msg=>"job_storage_dir: $job_storage_dir") ;
  }
  my $job_dir = $job_storage_dir.'/'.$job_id ;
  if (! -d $job_dir) {
    SrnaTools::Exception::FileAccess->throw(message=>"The job $job_id does not exist on the server. The URL might be incomplete or your job has already been deleted. Please note that we store job data for one week only.",log_msg=>"Missing job dir: $job_dir (may have expired)") ;
  }

} # _job_dir_health_check


#######################################################
# Get the number of jobs on the cluster that are
# currently in the queue (includes those that are waiting
# to be started).
#######################################################
sub _number_jobs_in_remote_queue{
  my $self = shift;
  my $n_jobs;
  
  my $user_host_name = $self->user_host_name('cluster') ;
  my ($user_name) =($user_host_name=~/^(.+?)@/);
  
  SrnaTools::Exception::Misconfiguration->throw("user_host_name string not defined in config or user name could not be parsed") unless $user_name;
  
  my $cmd = "ssh ".$user_host_name." qstat -u ".$user_name.' 2>/dev/null | wc -l 2>/dev/null' ;
  eval {
    $n_jobs = `$cmd` ; 
  };
  if ($@ || ! defined $n_jobs) { # no response
    SrnaTools::Exception::Backend::RemoteHostNoResponse->throw ;
  }
  
  $n_jobs = $n_jobs -2 if $n_jobs > 0 ; # remove the header lines 
  
  return $n_jobs;
}


#######################################################
# Check the downtime file to see whether we are in
# a server downtime. Return end of downtime if there 
# is one.
#######################################################
sub _within_scheduled_downtime{
  my $self = shift ;
  my $downtime_schedule = $self->root_dir('server').'/'.$self->cfg->{downtime_schedule_file};
  return unless $downtime_schedule && -r $downtime_schedule && !-z $downtime_schedule;
  
  use Time::Piece;
  my $time_now = localtime;
  my $max_end_time ;
  open (SF, '<', $downtime_schedule) || return ;
  while (<SF>) {
    chomp ;
    my ($begin_time,$end_time) ;
    /^(.+?) (.+)$/ ;
    eval{
      $begin_time = Time::Piece->strptime($1, "%d-%m-%Y;%H:%M");
      $end_time = Time::Piece->strptime($2, "%d-%m-%Y;%H:%M");
    } ;
    if ($begin_time <= $time_now && $end_time >= $time_now) {
      $max_end_time = $end_time if !$max_end_time || $end_time > $max_end_time ;
    }
  }
  close SF ;
  return $max_end_time ? $max_end_time->strftime("%d/%m/%y @ %H:%M (GMT)") : 0;
  
} # _within_scheduled_downtime


#######################################################
# Check if a job has already been sent to queue directory
#######################################################
sub _job_sent_to_queue_directory{
  my $self = shift ;
  my $job_id = shift ;
  
  if (! $job_id) {
    SrnaTools::Exception->throw("subroutine called with missing parameter job_id");
  }
  
  my $status_file = $self->path_to('job_storage_dir').'/'.$job_id.'/status/sent_to_queue_directory' ;
  return -e $status_file ? 1 : 0 ;
 
}


#######################################################
# Check if a job has been completed
#######################################################
sub _job_complete{
  my $self = shift ;
  my $job_id = shift ;
  
  if (! $job_id) {
    SrnaTools::Exception->throw("subroutine called with missing parameter job_id");
  }
  
  my $status_file = $self->path_to('job_storage_dir').'/'.$job_id.'/status/completed' ;
  return (-e $status_file) ? 1 : 0 ;
 
}


#######################################################
# Retrieve error msg (if any) from JOBDIR/status/errors
# file.
#######################################################
sub _job_error_msg{
  my $self = shift;
  my $job = $self->job ;
  SrnaTools::Exception->throw("Mo Job associated with App") unless $job;
  
  my $error_file = $self->path_to('job_storage_dir').'/'.$job->job_id.'/status/errors' ;
  if (-e $error_file){
    my $msg = `cat $error_file 2>/dev/null` ;
    $msg ||= "unknown exception - no details available";
    return $msg ;
  } else {
    return ;
  }
} #_job_error_msg


#######################################################
# Get the time of submission of this job
#######################################################
sub _get_job_submission_time{
  my $self = shift ;
  my $job_id = shift ;
  my $time = '...' ; # default
  
  if (! $job_id) {
    SrnaTools::Exception->throw("subroutine called with missing parameter job_id");
  }
   my $flag_file = $self->path_to('job_storage_dir').'/'.$job_id.'/status/sent_to_queue_directory' ;
  
  # Get the time of creation of this file
  if (-r $flag_file) {
    eval{
    my $stats = stat($flag_file);
    $time = scalar localtime $stats->mtime;
    } ;
  }
  return $time ;
} # _get_job_submission_time


#######################################################
# The cgi->Vars method returns a hash of parameters from
# the URL string but it squeezes multi-valued parameters
# into a single nullcharacter-separated string. Here we
# turn those into arrays. There doesn't seem
# to be a way of achieving this with CGI.pm directly.
# Returns a reference to the parameters hash.
# Also remove all potentially dangerous characters
# and untaint user input. 
#######################################################
sub _parse_url_params{
  my $self = shift ;
  my %vars = $self->query->Vars ;
  
  foreach my $key ( keys %vars ) {
    my $value = $vars{$key};
    if ($value=~m/\0/){ # multivalue (null-char separated)
      my @values ;
      foreach ( split( /\0/, $value) ) {
        push @values, $self->_sanitize_param($_) ;
      }
      $vars{$key} = \@values;
    } else { # just a single value
      $vars{$key} = $self->_sanitize_param($value) ; 
    }
  }
  
  return \%vars ;
} # _parse_url_params


#######################################################
# Remove all unsafe characters from user input. This
# includes any double-dots to prevent directory changes
# but single dots must be allowed for floating point numbers.
#######################################################
sub _sanitize_param{
  my ($self,$param_value) = @_ ;
  my $safe_characters = '>\n a-zA-Z0-9_\-\.@';
  
  return unless defined $param_value ;
  
  #$param_value =~s/\s/_/g ; # replace whitespace with underscore
  $param_value=~s/\.{2,}/./g; # replace multiple dots with one
        
  # remove all unsafe characters   
  $param_value=~ s/[^$safe_characters]//g;
  # now we can untaint it
  if ($param_value=~/^([$safe_characters]*)$/) {
    return $1 ;
  } else {
    SrnaTools::Exception->throw(message=>"$param_value contains only illegal characters. Please make sure all of your input consists of letters and numbers only.",log_msg=>"param_value: $param_value") ;
  }
} # _sanitize_param


#######################################################
# upload files to /data
# in job data sub-dir.
# Use File::Copy move (which basically uses the
# 'rename' command if possible) to move files from the
# cgi temp dir to the job directory.  
# Usage:
# $self->_upload_file( cgi_field_name, destination_dir)
#######################################################
sub _upload_file{
  my $self = shift ;
  my $field_name = shift ;
  my $dest_dir = shift ;
  my $q = $self->query;
  
  SrnaTools:Exception->throw(message=>"missing arguments: field_name and dest_dir required",log_msg=>"field_name: $field_name, dest_dir: $dest_dir") unless $field_name && $dest_dir;
  
  # Get the actual file name of the CGI-
  # generated tmp file and move it
  my $cgi_upload_param = $q->param($field_name) ;
  
  # just give up on ampty field - input 
  # validation must be done by the Job
  return unless $cgi_upload_param ; 
  
  my $tmp_filename = $q->tmpFileName($cgi_upload_param); 
  
  my ($untainted_filename) = ($tmp_filename=~/^(.+)$/);
  if (-z $untainted_filename ) {
    SrnaTools::Exception::Fileupload->throw(message=>"uploaded file was empty.",log_msg=>"tmp_filename: $tmp_filename") ;
  }
  
  my $target_file = $dest_dir.'/'.$field_name.'.raw' ;
  move( $untainted_filename, $target_file) || SrnaTools::Exception::Fileupload->throw(message=>"Could not upload file",log_msg=>"target_file: $target_file") ;
  
} # upload_file


#######################################################
# Get disk usage from cluster
# or local machine, depending
# on whether we are running in
# local batch or cluster mode.
# Run du -s (summary) through
# ssh if running on cluster
#######################################################
sub _get_job_dir_usage{
  my $self = shift ;
  my $execution_env = shift ;
  
  my $usage ;
  my $cmd ;
  my $output ;
  
  my $job_path = $self->path_to('queue_job_dir',$execution_env) ;
  if (!$job_path) {
    SrnaTools::Exception::Misconfiguration->throw("queue_job_dir for exectuion environemt '$execution_env' not defined in config") ;
  }
  
  if($execution_env eq 'server' or $execution_env eq 'local_queue') { # running locally
    $cmd = "du -s ".$job_path;
  } else { # running on cluster
    my $user_host_name = $self->user_host_name('cluster') ;
    if (!$user_host_name) {
      SrnaTools::Exception::Misconfiguration->throw("user_host_name string not defined in config") ;
    }
    $cmd = "ssh ".$user_host_name." du -s ".$job_path.' 2>/dev/null' ;
  }
  eval {
    $output = `$cmd` ; 
  };
  
  if ($@ || !$output) { # no response
    SrnaTools::Exception::Backend::RemoteHostNoResponse->throw ;
  }
  
  ($usage) = ($output=~/^(\d+)/) ; 
  return $usage ;

} # get_job_dir_usage


#######################################################
# Handle errors
# Parse the error name and details according to debug
# level. Write to error-log and forward to error
# run mode
#######################################################
sub handle_error{
  my ($self, $err) = @_ ;
  
  my ($err_name, $details, $log_details) = $self->parse_error($err) ;
  warn "### ERROR ### Caught a $err_name error. Details: $log_details" ;
  $self->forward('error_page', $err_name, $details) ;
}


#######################################################
# access web-specific config data
#######################################################
# direct access
sub title{ $_[0]->cfg->{title} || '' }
sub sidenav_partial{ $_[0]->cfg->{sidenav_partial} } # TODO misconfig error checking


#######################################################
# helpers used in templates
#######################################################
# link to a document with the name of the doc 
# doesn't need '.tt' suffix, the link text and the id
# on the page to scroll to (may contain the leading #)
# but not required
sub static_doc_link{
  my $doc = shift;
  my $text = shift || $doc;
  my $id = shift || '';
  $doc=~s/\.tt$// ; 
  $id = '#'.$id if $id && $id!~/^#/;
  return '<a href="'.THIS_SCRIPT_LINK.'?rm=document&page='.$doc.$id.'">'.$text.'</a>';
}
# same but shortcut for a link with a small question mark
sub question_mark_link{
  my ($doc, $id) = @_ ;
  my $text = '(?)' ; # could replace this with image
  my $html = static_doc_link($doc, $text, $id);
  $html=~s/">/" class="qmarkLink">/ ;
  return $html ;
}

# link to form page of a tool
sub tool_link{
  my $tool = shift;
  my $text = shift || $tool;
  return '<a href="'.THIS_SCRIPT_LINK.'?rm=input_form&tool='.$tool.'">'.$text.'</a>';
}

# this is just the src of the image, not a
# full image tag
sub image_src{
  my $file = shift || return '' ;
  return '../images/'.$file;
}
# link tag to an examples file
sub files_link{
  my $file = shift || return '' ;
  my $text = shift || '';
  return '<a href="../files/'.$file.'">'.$text.'</a>';
}


# DO NOT DELETE
1;
