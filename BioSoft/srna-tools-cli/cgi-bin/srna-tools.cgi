#!/usr/bin/perl -T

##########################################################################
# 
# This is the runner script for the web-based version of the UEA sRNA
# toolkit. It runs a CGI::Application.
#
##########################################################################

use strict ;
use warnings ;
use CGI qw/:standard/ ;
use CGI::Carp qw(fatalsToBrowser set_message) ;
#use FindBin qw($Bin); # find script directory
#use lib "$Bin/../lib"; # not required for cgi script
use lib '../lib'; # The SrnaTools are in the lib dir
use SrnaTools::Application::WebApp;
my $webapp = SrnaTools::Application::WebApp->new();

$webapp->run();

