#!/usr/bin/perl

#####################################################################
#
# This is the command line version of the UEA sRNA toolkit
#
#####################################################################

use strict;
use warnings;
use File::Copy ; 
#use Getopt::Param; # cgi-like param parsing
use FindBin qw($Bin); # find script directory
use lib "$Bin/lib"; 
use SrnaTools::Application::CliApp;


#######################################################
#  some Declarations
#######################################################
my $APP_CONFIG_FILE = "$Bin/config/application.conf";
my $TOOLS_CONFIG_FILE = "$Bin/config/tools.conf";

# Instantiate SrnaTools application
my $app = SrnaTools::Application::CliApp->new(
  cfg_file => $APP_CONFIG_FILE,
  tools_cfg_file => $TOOLS_CONFIG_FILE, 
) ;

$app->run ;



exit 0 ;
