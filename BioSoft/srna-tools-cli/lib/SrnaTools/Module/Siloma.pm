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
# It returns a gbrowse image, a fasta file of matching sequences
# and a stats file for sRNAs matching a backbone sequence
# 
# 
# TODO This class has been converted from a standalone script and might
# still need some cleanup (e.g. passing of parameters)
# OLD DOC
# Arguments
# --working_dir   Working-directory for this project
# --plot-nr       plot arrows in non redundant form with
#                 arrow thickness representing the log10 of the abundance
# --region_length length of region to plot. This can't be extracted from the 
#                 Patman output file
# --plot-labels   plot labels at arrows (sRNA sequence)

package SrnaTools::Module::Siloma ;
use base SrnaTools::Module ;
use strict ;
use warnings;
use File::Copy ;

sub run {
my $self = shift;
require Bio::Graphics;
require Bio::SeqFeature::Generic;


#############################
# Parameters and
# declarations
#############################
my $module_name = "SiLoMa" ; # for error log

my $working_dir = $self->job->job_working_dir;
my $result_dir ;
my $job_name = $self->job->job_name;
my $pat_out_file; 
my $match_seqs_file ;
my $data_dir ;
my $region_length = $self->job->cache->{region_length};
my $plot_nr = $self->param('plot_nr');
my $plot_labels = $self->param('plot_labels') ;
my $backbone_file;
my $max_num_arrows = $self->job->config->{max_num_arrows};

#############################
# Check parameters and files
#############################

eval {
  die "missing working directory parameter\n" unless $working_dir ;
  die "missing job name parameter\n" unless $job_name ;
  die "could not fetch region length from cache" unless $region_length;
  
  # Default error file is WORKING_DIR/errors
  $data_dir = $working_dir.'/data/' ;
  $backbone_file = $data_dir .'/'. $self->param('backbone_file') ;
  $pat_out_file = $data_dir.'/'. $self->param('pat_out_file') ;
  
  
  # Check files and directories
  die "could not find/read working directory\n" unless -d $working_dir ;
  die "data directory not found in working directory\n" unless -d $data_dir ;  
  die "could not find/read Patman result file\n" unless -r $pat_out_file ;
  
  die "could not find zip executable - please contact administrator\n" unless $self->binary_in_path('zip');
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error fetching/writing backbone sequence: $@") ;
}

#############################
# Initialize genome browser view
# with one track for the sRNA arrows
# and one for the ruler
#############################
# update status on server
$self->update_status("$module_name: initialising genome browser image") ;
my $panel ;
my $srna_track ;
eval {
  $panel = Bio::Graphics::Panel->new(
    -length => $region_length,
    -width  => 800,
    -pad_left => 10,
    -pad_right => 10,
  );
  die "failed to initialise Bio::Graphics::Panel" unless $panel;
  
  # Add ruler
  my $full_length = Bio::SeqFeature::Generic->new(-start=>1,-end=>$region_length);
  die "failed to add a sequence feature" unless $full_length;
  $panel->add_track(
    $full_length,
    -glyph   => 'arrow',
    -tick    => 2,
    -fgcolor => 'black',
    -double  => 1,
  );
  
  # Add track for sRNAs
  # The line width and height are fixed if we are plotting
  # in redundant form or dynamic (log10 scale) if plot-nr
  # option is used
  my $line_width = 2 ;
  my $line_height = 6 ;
  if ($plot_nr) {
    $line_width = \&line_width_from_score ;
    $line_height = \&line_height_from_score ;
  }
  my $labels = $plot_labels ? 1 : 0 ;
  $srna_track = $panel->add_track( 
    -glyph      => 'arrow',
    -label      => $labels ,
    -fgcolor    => \&fgcolor,
    -linewidth  => $line_width,
    -height     => $line_height,
  );

} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error setting up genome browser image: $@") ;
}

#############################
# parse patman file and add arrow to
# genome browser view for each match
# and also add matching sequences to 
# FASTA file
#############################
# update status on server
$self->update_status("$module_name: generating output files") ;
eval {
  
  my $dest_dir = $working_dir.'/results/mapping_results' ;
  mkdir $dest_dir or die "could not generate result directory\n" ;
  
  open(PAT_OUT, $pat_out_file) or die "could not open PatMan output file\n" ;
  
  my %hits ; # for hit stats
  my $arrow_count;
  while ( my $patline = <PAT_OUT>) {
    # the sRNA IDs are in the format "SEQUENCE_count"
    # All hits are to the same ref seq so we don't need to parse
    # its ID here
    my ($seq, $seq_count, $start, $stop, $strand) ;
    if ($patline=~/^.+?\t([AGCTU]+)_S[\d+]\((\d+)\)\t(\d+)\t(\d+)\t([\+-])/  ) {
      ($seq, $seq_count, $start, $stop, $strand) = ($1,$2,$3,$4,$5) ;
      my $size = length $seq ;
      $hits{nr}{seq_match_count}{$seq} = 1;
      $hits{nr}{seq_size_count}{$size}{$seq} = 1 ;
      $hits{nr}{seq_strand_count}{$strand}{$seq} = 1 ;
      $hits{nr}{sequence}{$seq}{count} = $seq_count ;
      push @{$hits{nr}{sequence}{$seq}{pos}}, $start.'..'.$stop.$strand ;
      $hits{red}{match_count}+= $seq_count ;
      $hits{red}{size_count}{$size}+= $seq_count ;
      $hits{red}{strand_count}{$strand}+= $seq_count ;
      
      if ($plot_nr){
        $arrow_count++ ;
      } else {
        $arrow_count += $seq_count;
      }
      
      if ($max_num_arrows && ($arrow_count > $max_num_arrows)){
        my $msg = "$module_name, exceeded maximum number of features to draw: tried to draw $arrow_count features but maximum is $max_num_arrows. Try to select a smaller region.";
        $msg .=" Currently trying to draw one arrow for each match: try the option 'Plot sRNA hits in non-redundant form' instead." unless $plot_nr;
        SrnaTools::Exception::ModuleExecution->throw($msg) ;
      }
      
    } else {
      die "patman result line $patline could not be parsed\n" ;
    }  
    
    # Add feature to track
    # uncomment the loop and replace linewidth and height with 
    # fixed values to print mutliple arrows for redundant seqs
    my $feature = Bio::SeqFeature::Generic->new(
      -display_name=> $seq.'_'.$seq_count,
      -start       => $start,
      -end         => $stop,
      -strand      => $strand,
      -score       => $seq_count,
    );
    die "failed to add feature to genome browser image" unless $feature;
    
    # if we are doing redundant format
    # we plot this arrow x times where x=seq_count
    my $plot_times = $plot_nr ? 1 : $seq_count ;
    foreach (1..$plot_times) {
      $srna_track->add_feature($feature);
    }
  
  } # parser (patman file)
  close PAT_OUT;

  # generate image
  my $image_file = $dest_dir.'/'.$job_name.'_image.png' ;
  open (IMAGE_FILE, '>', $image_file) || die "Could not generate image file\n" ;
  print IMAGE_FILE $panel->png;
  close IMAGE_FILE ;
  
  # print non-redundant matching seq to fasta
  my $match_seqs_file = $dest_dir.'/'.$job_name.'_matches.fasta';
  open (SEQS_FILE, '>', $match_seqs_file) || die "Could not generate matching sequence fasta file\n" ;
  my $i = 1 ;
  foreach my $seq (sort keys %{$hits{nr}{sequence}} ) {
    my $count = $hits{nr}{sequence}{$seq}{count} ;
    my $pos_str = join(';',@{$hits{nr}{sequence}{$seq}{pos}}) ;
    print SEQS_FILE '>'.$i.'_'.$count.'x_pos:'.$pos_str."\n".$seq."\n" ;
    ++$i ;
  }
  close SEQS_FILE ;
  
  # copy backbone file into result dir
  my $ref_seq_dest = $dest_dir.'/'.$job_name.'_reference-sequence.fasta' ;
  copy($backbone_file,$ref_seq_dest) or die "Failed to copy reference sequence file\n" ; 

  # write stats to file
  my $stats_file = $dest_dir.'/'.$job_name.'_stats.txt' ;

  open (STATS_FILE, '>', $stats_file) || die "Could not generate stats file\n" ;
  
  print STATS_FILE "\"Overview of matching sRNAs for $job_name\"\n" ;
  print STATS_FILE "\"Read counts\"\n";
  print STATS_FILE "\"\",\"total\",\"non-redundant sequence count\"\n" ;
  
  my $nr_match_count = keys %{$hits{nr}{seq_match_count}} ;
  my $red_match_count = $hits{red}{match_count} || 0 ;
  print STATS_FILE "\"All matches:\",$red_match_count,$nr_match_count\n" ;

  # strands 
  foreach my $strand (sort keys %{$hits{red}{strand_count}} ) {
    my $nr_count = keys %{$hits{nr}{seq_strand_count}{$strand}} ;
    my $red_count = $hits{red}{strand_count}{$strand} || 0 ;
    print STATS_FILE "\"Matches $strand strand\",$red_count,$nr_count\n" ;
  }
  
  # sizes
  foreach my $size (sort {$a <=> $b} keys %{$hits{red}{size_count}} ) {
    my $nr_count = keys %{$hits{nr}{seq_size_count}{$size}} ;
    my $red_count = $hits{red}{size_count}{$size} ;
    print STATS_FILE "\"Matches of length $size nt\",$red_count,$nr_count\n" ;
  }
  close STATS_FILE ;
  
  # zip results: -q quite, -r recursive, -j exclude full path to dir
  my $cmd = "zip -qrj $dest_dir $dest_dir" ;
  system($cmd) == 0 or die "Failed to compress result files\n" ; 
  system("rm -rf $dest_dir") ; # remove original directory
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while generating output files: $@") ;
}
return 1;
} #run


sub fgcolor {
  my $feature=shift;
  my $len = ($feature->end - $feature->start +1) ; 
  return "pink" if $len >=15 and $len <20; 
  return "red" if $len >= 20 and $len <=21; 
  return "green" if $len >=22 and $len <=23 ;
  return "blue" if $len >=24 and $len <=25 ; 
  return "gray" ;
}
  
sub line_width_from_score {
  my $feature=shift; 
  my $score = $feature->score; 
  if ($score){
    return ((log($score)/log(10))*2)+2 
  } else{
    return 2
  } 
}
sub line_height_from_score {
  my $feature=shift; 
  my $score = $feature->score; 
  if ($score){
    return ((log($score)/log(10))*6)+6 
  }else{
    return 6
  } 
}

1;