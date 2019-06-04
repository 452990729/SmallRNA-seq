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
# This class contains definitions for Exceptions that the SrnaTools
# modules can use. It uses the Exception::Class module
# Each exception also accepts a filed 'log_msg' which is only printed
# to the error log, except if the debug option is used in the 
# application config.
#
#####################################################################

package SrnaTools::Exceptions;
use strict;
use warnings;
use Exception::Class(
  'SrnaTools::Exception' => { 
    description => 'Generic base class for Srna Tools exceptions',
    fields => [ 'log_msg' ], 
  },
  'SrnaTools::Exception::Misconfiguration' => {
    isa         => 'SrnaTools::Exception',
    description => 'Missing item in configuration file',
  },
   'SrnaTools::Exception::FileAccess' => {
    isa         => 'SrnaTools::Exception',
    description => 'A file does not exist or is not read/writable',
  },
  'SrnaTools::Exception::ToolNotAvailable' => {
    isa         => 'SrnaTools::Exception',
    description => 'Throw when a tool was requested that is'.
                   ' not availvable as per site config file'
  },
  'SrnaTools::Exception::Downtime' => { 
    description => 'Throws when we are within server downtime',
    isa         => 'SrnaTools::Exception',
  },
  'SrnaTools::Exception::TemplateError' => {
     isa         => 'SrnaTools::Exception',
  },
   'SrnaTools::Exception::Url' => {
     isa         => 'SrnaTools::Exception',
     description => 'used for cases where a resource does not exist'.
                    'or URL parameter is missing'
  },
  'SrnaTools::Exception::Url::UrlParameterMissing' => {
     isa         => 'SrnaTools::Exception::Url',
  },
  'SrnaTools::Exception::Url::ResourceNotAvailable' => {
     isa         => 'SrnaTools::Exception::Url',
  },
  'SrnaTools::Exception::Backend' => {
    isa         => 'SrnaTools::Exception',
    description => 'Base for errors on the (remote) backend'
  },
  'SrnaTools::Exception::Backend::RemoteHostNoResponse' => {
    isa => 'SrnaTools::Exception::Backend'
  },
  'SrnaTools::Exception::Fileupload' => {
    isa => 'SrnaTools::Exception::Backend'
  },
  'SrnaTools::Exception::ExceedJobDirQuota' =>{
    isa => 'SrnaTools::Exception::Backend'
  },
  'SrnaTools::Exception::ExceedJobNumberQueue' =>{
    isa => 'SrnaTools::Exception::Backend'
  },
  'SrnaTools::Exception::JobPreparation' => {
    isa         => 'SrnaTools::Exception',
  },
  'SrnaTools::Exception::Parameter' => {
    isa         => 'SrnaTools::Exception',
  },
   'SrnaTools::Exception::ModuleExecution' => {
    isa         => 'SrnaTools::Exception',,
    description => 'For errors in the actual back-end scripts'
  },
);

1;