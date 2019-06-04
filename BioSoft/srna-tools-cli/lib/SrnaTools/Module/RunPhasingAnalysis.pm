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
# It runs the ta-siRNA analysis
# TODO This class has been converted from a standalone script and might
# still need some cleanup (e.g. passing of parameters)
#####################################################################

package SrnaTools::Module::RunPhasingAnalysis ;
use base SrnaTools::Module ;
use strict;
use warnings;



sub run{
my $self = shift ;

my $input_file = $self->param('pat_out_file');
my $genome_file = $self->param('genome');
my $pval = $self->param('pval');
my $min_abundance = $self->param('abundance');

# update status on server
$self->update_status("running phasing analysis") ;

if (! -r $input_file){
  SrnaTools::Exception::ModuleExecution->throw(message=>"Could not read genome matches file",log_msg=>"file:$input_file");
}
$self->job->_run_module(
  'phasing',{
    input_file  => $input_file,
    genome_file => $genome_file,
    pval        => $pval,
    min_abundance => $min_abundance,
  }
);

return 1;
}


1;