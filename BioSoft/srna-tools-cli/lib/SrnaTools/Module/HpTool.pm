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


####################################################################
#
# This is a module of the sRNA tool kit.
# It rus RNAfold on the long sequence, maps the short 
# sequences onto the long one and highlights them on the secondary
# structure plot.

package SrnaTools::Module::HpTool ;
use base SrnaTools::Module ;
use strict ;
use warnings;
use File::Temp qw( tempfile tempdir );
use IO::String ;
use Cwd;

sub run{
my $self = shift;

require Bio::SeqIO ;

# Commands 
my $RNAfold_cmd = 'RNAfold' ; 
my $RNAplot_cmd = 'RNAplot' ;
my $ps2pdf_cmd = 'ps2pdf' ; 
my $ps2eps_cmd = 'ps2eps' ;
my $working_dir =  $self->job->job_working_dir;
my $long_seq_param = $self->param('longSeq');
my $short_seqs_param = $self->param('shortSeqs');
my $cwd = getcwd;
($cwd) =($cwd=~/^(.+)$/);

# doc-root relative paths to temp dir
if (! -d $working_dir) {
  SrnaTools::Exception::ModuleExecution->throw(message=>"Could not access job working directory",log_msg=>"trying to access: $working_dir") ;
}

# rgb colors for highlighting short seqs
my @rgb_colors = ('255 40 40', 
                  '255 0 255', 
                  '0 0 255', 
                  '0 255 255', 
                  '0 255 0', 
                  '255 255 0', 
                  '255 127 0',
                  '255 184 235',
                  '181 255 184',
                  '199 255 176',
                  '127 62 62',
                  '255 227 99',
                  '127 0 0', 
                  '127 0 127', 
                  '0 0 127',
                  '0 127 127',
                  '0 127 0',
                  '130 127 0',
                  '153 153 153',
                  '204 204 204') ;

  # parse long sequence
  my %long_seqs = () ;
  parse_sequences( $long_seq_param, \%long_seqs) ;
  my ($long_seq_id) = keys %long_seqs ;
  my $number_long_seqs = keys %long_seqs ;
  my $long_seq = $long_seqs{$long_seq_id} ;
  if ($number_long_seqs > 1) {
    SrnaTools::Exception::ModuleExecution->throw("More than one long sequence to fold was entered. Please enter a single long sequence") ;
  } elsif ( !$number_long_seqs ) {
    SrnaTools::Exception::ModuleExecution->throw("No long sequence to fold or incorrect format (not valid FASTA)") ;
  }
  
  # parse short sequences
  my %short_seqs = ();
  parse_sequences( $short_seqs_param, \%short_seqs ) ;
  
  # It is ok to run this without any short sequences but if sequences
  # are given we need to check that they could be read
  if (! keys %short_seqs ) {
    SrnaTools::Exception::ModuleExecution->throw("Some or all of the short sequences are not in valid FASTA format") ;
  }
  
  ### get match positions of short sequences and
  # construct postscript macro for RNAplot
  # in the format 
  # "start_pos end_pos 8 color omark" (8 is the thickness of the line)
  # start the macro by putting a circle around the
  # 5' end
  my $rnaplot_ps_macro = '1 cmark ' ;
  
  my $i = 0 ;
  foreach my $short_seq_id (sort keys %short_seqs){
    my @positions = get_pos(\$long_seq, \$short_seqs{$short_seq_id}) ;
    foreach (@positions) {
      $rnaplot_ps_macro .= $_ .' 8 ' . convert_rgb_to_ps_format($rgb_colors[$i]) . ' omark ';
    }
    ++$i;
    last if $i > 19 ;
  }
  
  ### run RNAfold
  my $RNAfold_out_file = $working_dir.'/data/RNAfold_out';

  my $cmd = "echo '$long_seq' | $RNAfold_cmd > $RNAfold_out_file" ; # run RNAfold
  my $out = `$cmd` ;
  # There should not be an output of this command, if so we have an error
  if ($out) {
    SrnaTools::Exception::ModuleExecution->throw("There was a problem running RNAfold") ;
  }

  ### run RNAplot
  # RNAplot always generates its output in a file named after the input sequence
  # or (if not fasta) "rna.ps". Therefore we have to generate another directory
  my $RNAplot_out_file = 'rna.ps' ; # default RNAplot file name
  my $RNAplot_ps_final = 'rna.ps_final' ; # file after addition of label for pdf 

  # make temp dir
  my $RNAplot_dir = $working_dir.'/data/RNAplot_dir';
  
  mkdir $RNAplot_dir;
  chmod 0777, $RNAplot_dir ; 
  chdir $RNAplot_dir ; 

  # RNAplot command
  $cmd = "$RNAplot_cmd --pre \"$rnaplot_ps_macro\" < $RNAfold_out_file " ; 
  $out = `$cmd` ;
  # There should not be an output of this command, if so we have an error
  if ($out) {
    SrnaTools::Exception::ModuleExecution->throw("There was a problem running RNAplot") ;
  }

  ### modify and convert images
  # We add a label to the postscript file, then convert it to pdf
  # and also generate the jpg file to print to the result page
  
  # add label to postscript file
  # open the RNAplot outputfile in its temp dir and another file for writing the modified ps
  my $label ="Secondary structure for '$long_seq_id'" ;

  open PS , '<', $RNAplot_dir .'/'. $RNAplot_out_file or SrnaTools::Exception::ModuleExecution->throw("could not open RNAplot out file") ;
  open PS_FINAL, '>' , $RNAplot_dir . '/' . $RNAplot_ps_final or  SrnaTools::Exception::ModuleExecution->throw("could not generate result files") ;

  # add label to ps file
  # this is a bit of a hack:
  # the ps code is not inserted into the ideal position in the file
  # but by inserting just before "%%BeginProlog" we can print the label before the coordinate system is translated for plotting the RNA
  # post-script stack used:
  # Helvetica findfont -> change font
  # 12 scalefont -> set font size
  # 100 10 moveto -> set cursor position (bottom of page)
  # setfont
  # $label show -> print the label
  while (<PS>) {
    s/%%BeginProlog/\/Helvetica findfont\n12 scalefont\n100 10 moveto\nsetfont\n\($label\) show\n%%BeginProlog/ ;
    print PS_FINAL;
  }
  close PS ;
  close PS_FINAL ;

  # convert to PDF and write to results file
  my $pdf_file = $working_dir.'/results/Structure_plot.pdf';

  $cmd = $ps2pdf_cmd . ' '. $RNAplot_dir . '/' . $RNAplot_ps_final . ' ' . $pdf_file ;
  chmod 0666, $pdf_file ; 
  my $output = `$cmd` ;
  if ($output) {
    SrnaTools::Exception::ModuleExecution->throw("There was a problem running ps2pdf: $output") ;
  }
  # Make bitmap for displaying on page
  my $jpg_file = $working_dir.'/results/Structure_plot_bitmap.jpg';
      
  # generate jpeg file from ps using ghostscript
  # -r -> resolution
  # -dBATCH -dNOPAUSE -q -> silent, non-interactive mode
  # -dEPSCrop -> crop image to boundary box created by ps2eps
  # - (at the end) -> get input file from STDIN 
  $cmd = 'gs -sDEVICE=jpeg -r100x100 -dGraphicsAlphaBits=4 -sOutputFile='.$jpg_file.' -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -q - < '.  $RNAplot_ps_final ;
  $output = `$cmd` ;
  chmod 0666, $jpg_file ; 
  
  # create data for label and legend in results page
  my $legend = { label => $label};
  my $ii=0 ;
  foreach my $short_seq_id (sort keys %short_seqs){
    $legend->{short_seqs}{$short_seq_id} = get_css_rgb($rgb_colors[$ii]);
    ++ $ii;
    last if $ii > 19 ;
  }
  $self->job->_write_hash_to_file($legend,$working_dir.'/results/legend.txt');
  chdir $cwd; # back to where we came from
  return 1;
}

# Convert a space delimited rgb value into css notation and return a string
# rgb(R,G,B)
sub get_css_rgb {
  my $input_rgb = shift ; # in format R G B
  my ($r,$g,$b)=($input_rgb=~/(\d+)/g) ;
  return "rgb($r,$g,$b)" ;
}

# Convert space delimited RGB into notation for RNAplot macros
# where 255 = 1
sub convert_rgb_to_ps_format {
  my $input_rgb = shift ; # in format R G B
  my ($r,$g,$b)=($input_rgb=~/(\d+)/g) ;
  $r = $r / 255 ;
  $g = $g / 255 ;
  $b = $b / 255 ;
  return "$r $g $b" ;
}

# Extract all sequences in fasta format from input string
# arguments: input_string, ref-to-hash (for results), ref-to-error-array
sub parse_sequences{
  my $input = shift ;
  my $seqs_ref = shift ;
  
  # convert to filehandle and pass to SeqIO
#  my $string_2_fh = new IO::String($input) ;
eval {
	die " could not find/read input file $input \n" unless -r $input ;
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("error in parameters: $@") ;
}

  my $parser = Bio::SeqIO->new( -file => $input ,
                                -format => 'fasta') ;

  # parse sequences, add error if an ID was found twice
  # or a sequence is longer than 1kb
  while (my $seq_obj =$parser->next_seq) {
    my $id =$seq_obj->id ;
    if (exists $$seqs_ref{$id}) {
      SrnaTools::Exception::ModuleExecution->throw("The ID '$id' was given to more than one sequence. Please ensure that all ID lines are unique") ;
    }
        
    my $seq = $seq_obj->seq;
    
    if (length($seq) > 1000) {
      SrnaTools::Exception::ModuleExecution->throw("Sequence '$id' is longer than the maximum of 1kb");
    }
    
    if (!$seq) {
      SrnaTools::Exception::ModuleExecution->throw("Sequence '$id' did not have a recognisable sequence line. Please make sure that all sequences are entered in valid FASTA format");
    }
    
    if ($seq!~/^[AGCTUN]+$/i) {
      SrnaTools::Exception::ModuleExecution->throw("Sequence '$id' contained characters other than AGCTU or N.");
    }
    
    # Convert T to U
    $seq =~tr/Tt/Uu/;
    
    $$seqs_ref{$id} = uc $seq;
  }
  
}

# get match positions of short sequence on long sequence
sub get_pos{
  my ($long_seq_ref, $short_seq_ref) = @_ ;
  my @pos = () ;
  my $long_seq_len = length($$long_seq_ref) ;
  my $short_seq_len = length($$short_seq_ref) ;
  for (my $i = 1 ; $i <= $long_seq_len - $short_seq_len + 1 ; ++$i) {
    if (substr($$long_seq_ref, $i - 1, $short_seq_len) eq $$short_seq_ref) { # we have a match
    my $start = $i ;
    my $stop = $i + $short_seq_len - 1 ;
    push @pos, $start . ' ' . $stop ;
    }
  }
  return @pos ;
}

sub file_basename{
  my $path = shift ;
  return $path=~/.+\/(.+)$/;
}


1;
