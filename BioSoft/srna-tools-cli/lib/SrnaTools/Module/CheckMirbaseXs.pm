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
# This is a module of the sRNA tool kit
# It checks modified mirBASE files for the presence of added
# "XX" at both ends of each sequence and terminates if these
# are not found. This can happen after an update if this
# was not done using the update_mirbase.pl script, which
# automatically generates these files.
# Only the first sequence is checked
######################################################################

package SrnaTools::Module::CheckMirbaseXs ;
use base SrnaTools::Module ;
use strict;
use warnings;



sub run{
my $self = shift ;

# CPAN module might be in non-standard location
# see Module.pm for details of how we add to @INC
require Bio::SeqIO;

#############################
# Declarations
#############################
my $module_name = "Checking miRBASE files" ; # for error log
my $public_data_dir = $self->job->app->path_to('data_dir',$self->job->execution_env);  ;
my $working_dir = $self->job->job_working_dir;
my $check_file = $self->param('file');


#############################
# Check parameters and files
#############################
$self->update_status("checking miRBase files") ;
eval {
	die "missing working directory parameter\n" unless $working_dir ;
	die "missing public data directory parameter\n" unless $public_data_dir ;
  
  
	die " could not find/read working directory\n" unless -d $working_dir ;
	die " could not find/read public data directory\n" unless -d $public_data_dir ;
 die "no miRBASE file given\n" unless $check_file; 
} ;
if ($@) {
	SrnaTools::Exception::ModuleExecution->throw("$module_name, error in parameters: $@") ;
}

# Do the check on first sequence only
eval {
	my $full_path = $public_data_dir .'/'.$check_file ;
	die "miRBASE file not readable\n" unless -r $full_path ;
  my $in_seqio = Bio::SeqIO->new( -format => 'fasta', -file => $full_path ) || die " problem creating Bio::SeqIO object\n";
  
  my $inseq = $in_seqio->next_seq ;
  my $seq_str = $inseq->seq;
  if ($seq_str !~/^XX.+XX$/) {
		die "miRBASE file has not been prepared for overhanging matches\n" ;
  }
} ;
if ($@) {
	SrnaTools::Exception::ModuleExecution->throw("$module_name, error while checking mirBase file: $@") ;
}

return 1 ;
}#run

1; 