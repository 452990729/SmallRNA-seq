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
# It runs the mircat pipeline.
# TODO This class has been converted from a standalone script and might
# still need some cleanup (e.g. passing of parameters)
#####################################################################

package SrnaTools::Module::Mircat ;
use base SrnaTools::Module ;
use strict;
use warnings;
use Cwd;

sub run{
my $self = shift ;

my $module_name = "miRCat" ;
my $data_dir = $self->job->app->path_to('data_dir',$self->job->execution_env) ;

# The genome file is usually in the public data dir but
# if we are running in interactive mode on the CLI App then
# it is a complete path. The parameter genome_by_path is set
# in that case
my $genome_file ;
if ($self->job->param('genome_by_path')){
  $genome_file = $self->param('genome');
} else {
  $genome_file = $data_dir.'/'.$self->param('genome');
}


my $lib_dir = $self->job->app->path_to('lib_dir',$self->job->execution_env); 
my $working_dir = $self->job->job_working_dir ;
my $gff = $self->param('gff');
my $abundance = $self->param('min_abundance');
my $minsize = $self->param('minsize');
my $maxsize = $self->param('maxsize');
my $genomehits = $self->param('genomehits');
my $window_length = $self->param('window_length');
my $max_unique_hits = $self->param('max_unique_hits');
my $max_percent_unpaired = $self->param('max_percent_unpaired');
my $max_overlap_length = $self->param('max_overlap_length');
my $percent_orientation = $self->param('percent_orientation');
my $min_paired = $self->param('min_paired') ;
my $max_gaps = $self->param('max_gaps') ;
my $min_gc = $self->param('min_gc') ;
my $min_hairpin_len = $self->param('min_hairpin_len') ;
my $cwd = getcwd;

# update status on server
$self->update_status("$module_name: running miRCat pipeline (this might take a while)") ;


if (!$minsize){
  $minsize = 20;
}
if (!$maxsize){
  $minsize = 22;
}
if (!$genomehits){
  $genomehits = 16;
}

if ((!$genome_file) || (!$lib_dir) || (!$working_dir) || (!$abundance) || (!$minsize) || (!$maxsize) || (!$genomehits)){
   SrnaTools::Exception::ModuleExecution->throw(message=>"$module_name, insufficient input parameters",log_msg=>"got these: genome_file: $genome_file, lib_dir: $lib_dir, working_dir: $working_dir, abundance $abundance, minsize: $minsize, maxsize: $maxsize, genomehits: $genomehits");
}

chomp ($genome_file, $lib_dir, $working_dir, $gff, $abundance, $minsize, $maxsize, $genomehits);
  if (system("mkdir $working_dir/results/results/")!=0) { SrnaTools::Exception::ModuleExecution->throw(message=>"Could not make results directory",log_msg=>"tried to make:$working_dir/results/results/ ");
}

$self->job->_run_module(
  'process_hits',{
    genome_file          => $genome_file, 
    gff                  => $gff,
    min_abundance        => $abundance,
    minsize              => $minsize,
    maxsize              => $maxsize,
    genomehits           => $genomehits,
    window_length        => $window_length,
    max_unique_hits      => $max_unique_hits,
    max_percent_unpaired => $max_percent_unpaired,
    max_overlap_length   => $max_overlap_length,
    percent_orientation  => $percent_orientation,
    min_paired           => $min_paired,
    max_gaps             => $max_gaps,
    min_gc               => $min_gc,
    min_hairpin_len      => $min_hairpin_len,
    pval                 => $self->param('pval'),
    no_complex_loops     => $self->param('no_complex_loops'),
    hit_dist             => $self->param('hit_dist'),
    min_energy           => $self->param('min_energy'),
  } 
) ;
  
$self->update_status("$module_name: generating images") ;

$self->job->_run_module(
  'markup_hairpins',{
    input_file  => '/results/results/miRNA_hairpins.txt',
   }
);
eval{
  chdir "$working_dir/results/";
  system("zip -r -q results.zip results/ 1>/dev/null 2>&1")==0 or die;
};
if($@){
  SrnaTools::Exception::ModuleExecution->throw(message=>"$module_name, could not generate .zip archive from result files",log_msg=>"got this error: $@. Tried to cd into $working_dir/results/ and run'zip -r -q results.zip results/'");
}
system "rm -fr $working_dir/results/results/";
chdir $cwd ; # get out of here before we delete the dir (causes error)
return 1;

} ;

1;