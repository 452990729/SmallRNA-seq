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
# This is a module of the sRNA tool kit.

package SrnaTools::Module::FirepatCreateFailureOutfiles;
use base SrnaTools::Module ;

#use strict;
#use warnings;

sub run{
my $self = shift ;

my $moduleName = "FiRePat - create failure outfiles";
$output     = $self->job->job_working_dir.'/results/' ;
my $message = $self->param('message') || '-- sorry, could not identify problem --';

open out, '>',$output.'/no_clusters.txt' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output file");

print "FiRePat ran without errors but no cluster report could be generated from your data for the following reason:\n$message\n";

close out;

return 1;
}

1;