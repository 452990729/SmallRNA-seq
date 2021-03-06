######################################################################
# The UEA sRNA Toolkit: Perl Implementation - A collection of
# open-source stand-alone tools for analysing small RNA sequences.
# Copyright � 2011 University of East Anglia
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
# It fetches a requested region from one of the genome files
# and writes the sequence into /data/backbone_seq.fasta for
# the SiLoMa tool. 
# 
# TODO This class has been converted from a standalone script and might
# still need some cleanup (e.g. passing of parameters)
# OLD DOC
# Arguments
# --working_dir   Working-directory for this project
# --region        region in the genome in format CHROM:START..END
# --refseq        The path to the selected genome
# --genome_name   Display name of genome, not file name

package SrnaTools::Module::GetBackboneSeq ;
use base SrnaTools::Module ;
use strict ;
use warnings;


sub run{
my $self = shift;

require Bio::DB::Fasta ;
require Bio::SeqIO ;

#############################
# Parameters and
# declarations
#############################
my $module_name = "Fetch genomic region" ; # for error log

my $working_dir = $self->job->job_working_dir ;
my $public_data_dir = $self->job->app->path_to('data_dir',$self->job->execution_env);

# The genome file is usually in the public data dir but
# if we are running in interactive mode on the CLI App then
# it is a complete path. The parameter genome_by_path is set
# in that case
my $refseq ;
if ($self->job->param('genome_by_path')){
  $refseq = $self->param('refseq');
} else {
  $refseq = $public_data_dir.'/'.$self->param('refseq');
}


my $region_chrom = $self->param('region_chrom') ;
my $region_start = $self->param('region_start') ;
my $region_end = $self->param('region_end') ;
my $data_dir ;
my $output_file = $self->param('outfile') ;
my $genome_name = $self->param('genome_name');
my $pasted_seq_ref = $self->param('pasted_seq_ref') ;
my $region_length;

#############################
# Check parameters and files
#############################

eval {
  die "missing working directory parameter\n" unless $working_dir ;
  die "missing refseq parameter\n" unless $refseq ;
  if (! $pasted_seq_ref) {
    die "missing region parameter\n" unless $region_chrom && $region_start && $region_end;
    die "could not find/read genome file $refseq\n" unless -r $refseq ;
  }
  
  # Default error file is WORKING_DIR/errors
  $data_dir = $working_dir.'/data/' ;
  
  # Check files and directories
  die "could not find/read working directory\n" unless -d $working_dir ;
  die "data directory not found in working directory\n" unless -d $data_dir ;  
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error in parameters: $@") ;
}

# update status on server
$self->update_status("$module_name: fetching/writing backbone sequence") ;
    
eval {
  open (OUT, ">", $data_dir.$output_file) or die "Could not write to output file\n" ;
  
  if ($region_chrom ) {
    #############################
    # Fetch region from genome
    #############################
    my $db = Bio::DB::Fasta->new($refseq);
    my $seq_obj = $db->get_Seq_by_id($region_chrom);
    die "Could not find chromsome (BAC/scaffold) with ID $region_chrom in genome $genome_name. Please run the SiLoCo tool first to find the genomic regions of origin of your sRNAs.\n" unless $seq_obj ;
    
    my $len = $seq_obj->length ;
    die "The requested region start is outside chromosome $region_chrom for genome $genome_name (length:$len bp)\n" if $region_start > $len ;
    die "The requested region end is outside chromosome $region_chrom for genome $genome_name (length:$len bp)\n" if $region_end > $len ;
    
    my $subseq  = $seq_obj->subseq($region_start => $region_end);
    
    die "Could not extract region from chromosome $region_chrom for genome $genome_name\n" unless $subseq ;
    
    print OUT ">$region_chrom:$region_start..$region_end ($genome_name)\n$subseq\n" ;
    
    $region_length = $region_end - $region_start + 1;
 
  } else {
    #############################
    # write pasted seq to file
    #############################
    my ($id) = ($$pasted_seq_ref=~/^>(\S+)/) ? $1 : "pasted_sequence" ;
    my ($seq) = ( $$pasted_seq_ref=~/(\n[AGCTUN\s\n]+)$/i);
    $seq =~ s/[\n\s]//g ;
    $seq =~ tr/uU/tT/ ;
    $region_length = length ($seq) ;
    print OUT '>'.$id."\n".$seq ;
  }
}; #eval
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error fetching/writing backbone sequence: $@") ;
}
# we need the region length later
$self->job->cache->{region_length} = $region_length;
close OUT ;
return 1;
}

1;