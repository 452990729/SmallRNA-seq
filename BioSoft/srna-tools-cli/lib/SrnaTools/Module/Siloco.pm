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
# It runs the job of finding loci and
# comparing them.
#
# NOTE: this has been converted into a class from a non-OO script
# and still requires a bit of cleanup
#
# OLD DOCUMENTATION TODO needs to be updated
# Accept two input datasets (one of them can be a GEO sample)
# These files are preprocessed by uncompress_and_filter.pl
# and run through patman in run_patman.pl.
# The resulting Patman file contains sequence and abundance
# information for each of the two samples in the follwoing
# format in each Patman result line (seq ID field):
# SEQUENCE_NAME1(ABUNDANCE)[_NAME2(ABUNDANCE)]
# if a sequence occurs in only one of the samples, only
# that sample will be listed (i.e. abundance = 0 is not
# explicitly listed)
# The patman result file is expected to be
# WORKING_DIR/data/patman.out
#
# Input data:
# Two fasta sRNA files in non-redundant format
# generated by uncompress_and_filter.pl
# 
# Arguments
# --working_dir   Working-directory for this project
# -g              max gap length, defaults to 300bp
# -h              min weighted number of hits for a locus, defaults to 3
# --avg_rank_weight For weighting the rank sum (rank avg + rank ratio)
# --ratio_rank_weight 
# --links         Include Hyperlinks in csv file
# --sort_method   (pos|rank) Sort output by genome position (default) or
#                  by rank sum (lowest first = highest regulated/expressed)
# --pseucocount   Add this pseudocount to normalised read counts
# --error_file    This is where we log any errors from this script.
#                 It is a file ni the working dir that can be checked
#                 by the user later to see what went wrong.
#                 Default: WORKING_DIR/errors
# --single-sample Use a single sample and only get loci without comparing
# --show-unique    add a column for number of unqiquely matching reads
#
# The working directory is created by the cgi script that accepts the job
# it originaly contains the uploaded sample files in a subdir "data", which
# are mapped to the genome using Patman and are replaced by the result file
# called patman.out.
#
# We generate the output in the subdir "results" in the working directory. 
# Files:
# SAMPLE1_vs_SAMPLE2.csv      Contains all loci and their normalised weighted counts + log2 ratio
#                             and log10 avg abundancd
#
# Normalised count:
# First weight the number of reads by the number
# of matches to the genome, i.e. 100 reads at this
# position which also match to 100 other positions
# give a count of 1.
# Then normalise to the number of genome matching reads
# in the sample multiplied by a factor (1000)
# The final score is the -repetition-weighted count
# per 1 million genome-matching reads for a given sample

package SrnaTools::Module::Siloco ;
use base SrnaTools::Module ;
use strict ;
use warnings;
#use File::Basename ;

sub run{
my $self = shift ;

#############################
# Parameters and
# declarations
#############################
my $module_name = "SiLoCo" ; # for error log

my %srnas ; # nr list of input seqs
my %hits ; # stores the result of mapping sRNAs ro the genome
my %loci ; # stores the final loci
my @locus_ranks; # List of references into the %loci hash, sorted by rank sum

my $max_gap = $self->param('max_gap') ;
my $min_hits = $self->param('min_hits') ;

# For the ranking sum: the average and ratio
# ranks can be weighted seperately, e.g.
# give more emphasis on regulation than on 
# expression. 
# TODO: offer a 0.75 / 0.25 weighting in form
my $ratio_rank_weight = $self->param('ratio_rank_weight') ;
my $avg_rank_weight = $self->param('avg_rank_weight');

my @gb_links = @{$self->param('gb_links')};
my $working_dir  = $self->job->job_working_dir ;
my $data_dir ;
my @sample_names ;
my %samples ; # hash of sample names 
my $pseudocount = $self->param('pseudocount') ; # to avoid division by zero 
my %output_options ;
my $sort_method = $self->param('sort_method') ;
my $pat_out_file;
my $use_single_sample = $self->param('single_sample') ;
my $show_unique = $self->param('show_unique');


#############################
# Check parameters and files
#############################

# Check parameter sort method: default to "pos" if
# it isn't rank
$sort_method = 'pos' unless $sort_method eq 'rank' ;

eval {
  die " missing working directory parameter\n" unless $working_dir ;

  $data_dir = $working_dir.'/data/' ;
  $pat_out_file = $data_dir.$self->param('patman_out_file');
  
  # Check files and directories
  die " could not find/read working directory\n" unless -d $working_dir ;
  die " data directory not found in working directory\n" unless -d $data_dir ;  
  die "could not find/read Patman result file\n" unless -r $pat_out_file ;
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error in parameters: $@") ;
}

# Collect options for output
$output_options{gb_links} = \@gb_links ;
$output_options{sort_method} = $sort_method ;
$output_options{use_single_sample} = $use_single_sample ;
$output_options{show_unique} = $show_unique ;

# update status on server
$self->update_status("$module_name: started locus comparison") ;

#############################
# Parse hit positions from
# patman and build a hash
# of hit positions that
# is then expdanded into
# loci in the next step.
#############################

eval{ # Catch errors for this part
  $self->parse_hits( 
    \%srnas,
    \%hits, 
    $pat_out_file) ;
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while initialising loci: $@") ;
}

#############################
# Get the total number of genome
# matching reads for the samples
# for normalisation.
# Store in  
# $samples{$sample}{total_matching_reads}
#############################
$self->get_total_matching_read_counts(\%srnas, \%samples) ;

#############################
# Process hit maps and
# generate loci
# Combine hits based
# on proximity and
# Keep only loci that
# fulfill criteria.
#############################
# update status on server
$self->update_status("$module_name: building loci") ;

eval{ # Catch errors for this part
  $self->generate_loci(
    \%hits, 
    \%loci,
    \%samples, 
    $min_hits, 
    $max_gap
  ) ;
  
  # the hits and snras hashes aren't needed anymore
  %hits = () ;
  %srnas = () ;

} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while processing loci: $@") ;
}

#############################
# Get ratios and average expression
# for all loci to enable comparing
# them
# Skip this if single-sample
# is requested
#############################
# update status on server
$self->update_status("$module_name: comparing loci") ;
eval{ # Catch errors for this part
  $self->compare_loci(
    \%loci,
    \%samples,
    $pseudocount
  ) unless $output_options{use_single_sample};

} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while getting ratios and average expression levels: $@") ;
}


#############################
# Rank loci
#
# After "process_loci" we now have a list
# of above-threshold loci with their associated
# read counts for both samples.
# To give the user the best candidates we give each locus
# two ranks:
# One rank is based on the log2 ratios (equal ratios
# share a rank). The other rank is based on the average
# rexpression (normalised count)
# Then we can get a rank sum (rank A + rank B)
# an sort the results to find those loci that are both
# highly regulated and highly expressed..
# We can weight the two ranks to give e.g. more 
# weight to expression level than to regulation
# Skip if in single sample mode 
#############################
eval{ # Catch errors for this part
  $self->rank_loci(
    \%loci,
    \@locus_ranks,
    \%samples,
    $ratio_rank_weight,
    $avg_rank_weight,
  ) unless $output_options{use_single_sample};

} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while ranking loci: $@") ;
}

#############################
# Print results
#############################
# update status on server
$self->update_status("$module_name: generating result files") ;
eval{ # Catch errors for this part
  $self->generate_result_files(
    \%loci,
    $working_dir,
    \%samples,
    \%output_options,
  ) ;

} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while generating result files: $@") ;
}

return 1;
} # run


#############################
# parse_hit
# Parse the output from patman
# where abundances of each match
# sequence were stored in the
# ID line in the format
# SEQUENCE_sample1name(abundance)[_sample2name(abundance)
# Populate srna data with abundances
# and build a hash of hit positions
# that references the abundance 
# data.
#############################
sub parse_hits{
  my $self = shift ;
  my $srnas_ref = shift ;
  my $hits_ref = shift ;
  my $pat_out_file = shift ; # patman temp output file
	
  # Open the patman output file
  open(PAT, $pat_out_file) or die " could not open PatMaN output file\n" ;
  while ( my $patline = <PAT>) {
		
    # parse the hit
    # patman stores the entire ID line from
    # FASTA, not just ">\S+" like Bioperl
    my ($chrom_id, $seq, $seq_abund, $start, $stop, $strand) ;
    if ($patline=~/^(\S+).*?\t([AGCT]+)(_\S+)\t(\d+)\t(\d+)\t([\+-])/  ) {
      ($chrom_id, $seq, $seq_abund, $start, $stop, $strand) = ($1,$2,$3,$4,$5,$6) ;
    } else {
      die " Patman result line $patline could not be parsed\n" ;
    }
    # The abundances information is redundant, if a sequence
    # has more than one hit, each one will contain this
    # information. Check if we have the abundance information
    # for this sequence already and if not, add it by reading
    # into a hash of the sample names => abundance and
    # referencing to it from srnas hash.
    unless (defined $$srnas_ref{$seq}) {
      my %sample_abundance = ($seq_abund=~/_(\S+?)\((\d+)\)/g) ;
      my $n_samples = keys %sample_abundance ;
      die "More than two samples found\n" if $n_samples > 2;
      die "Sequence without sample data found\n" if $n_samples == 0 ;
      $$srnas_ref{$seq}{read_counts} = \%sample_abundance ;
    }
   
    # The values of the srnas hash are hashes with match and
    # read counts like this:      
    # {match_counts => x,
    # {read_count = { sample1 => x, sample2 => y}}
    # Here we link these to hit positions so that we can build
    # loci in the next step  
    my $seq_data = $$srnas_ref{$seq} ; 
    my $match_pos = $chrom_id.';'.$start.';'.($stop).';'.$strand ; # position in format "chrom_id;start;stop;strand"
    $seq_data->{match_count} += 1;
    $$hits_ref{$match_pos} = $seq_data ; # link the hit to the sequence data 
  } # parser (patman file)
  close PAT ;     
} # parse_hits
	
	
#############################
# get_total_matching_read_counts
# For normalising we need (for
# each sample) the count
# of reads with any match
# to the genome.
# Traverse all srnas hash
# which contains a match count
# and a read count for each
# non-redundant sequence.
# Add up read counts for all
# sRNAs that do have a match count
# To get the total number of reads
# matching the genome
sub get_total_matching_read_counts{
  my $self = shift ;
  my $srnas_ref = shift ;
  my $samples_ref = shift;
  
  foreach my $seq (keys %$srnas_ref) {
    next unless $$srnas_ref{$seq}{match_count} ; 
    foreach my $sample (keys %{ $$srnas_ref{$seq}{read_counts} } ) {
      $$samples_ref{$sample}{total_matching_reads} += $$srnas_ref{$seq}{read_counts}{$sample};
      $$samples_ref{$sample}{nr_matching_reads}++ ;
    } # foreach sample
  } # foreach seq
} # total_matching_read_counts



#############################
# generate_loci
# Traverse list of hits to the
# genome by position (lowest
# to highest start pos on each
# chromosome).
# Add hits to a locus if they
# are within the max gap.
#############################
sub generate_loci{
  my $self = shift ;
  my $hits_ref = shift ; # take hits from here
  my $loci_ref = shift ; # populate this hash
  my $samples_ref = shift ;
  my $min_hits = shift;
  my $max_gap = shift ;
  
  # Sort hits by chromosome and inside chromosome by start pos
  # Kepp track of the last chromosome_id and reset everything 
  # when starting on a new chromosome (loci can not span chromosomes)
  my $current_locus_start ; 
  my $current_locus_stop ;
  my $locus_num = 1 ; # counter for loci
  $$loci_ref{$locus_num} = {} ; # initialise first locus
  my $current_chrom = '' ;
  foreach my $hit (sort sort_by_chrom_and_start keys %{$hits_ref}) {
    # Extract positions
    my ($chrom,$start,$stop,$strand) = ($hit=~/^(.+?);(\d+);(\d+);([+-])$/) ;
    die "Could not read match position $hit\n" unless $chrom && $start && $stop && $strand ;
    
    # If the gap is too wide:
    # finalise the last locus and
    # keep if weighted total count 
    # above threshold. If not, delete
    # and re-use the locus number
    if (($chrom ne $current_chrom) || ($start > $current_locus_stop + $max_gap + 1) ) {
      # Keep the previous locus (locus num still points at that) ?
      if ($$loci_ref{$locus_num}{total_weighted_count} && $$loci_ref{$locus_num}{total_weighted_count} >= $min_hits ) { # keep it
        my $locus_pos = $current_chrom.';'.$current_locus_start .';'.$current_locus_stop ;
        $$loci_ref{$locus_num}{position} = $locus_pos ;
        ++$locus_num ;
        # must initialise new locus ref (otherwise add_to_locus below would fail)
        $$loci_ref{$locus_num}  = {};
      } else { # not above threshold - reset this locus
        $$loci_ref{$locus_num}  = {};
      } # if hits > threshold
      
      # Start the new locus with the current hit
      $current_locus_start = $start ;
      $current_locus_stop = $stop ;
      $current_chrom = $chrom ;
          
    } else { # still within gap limit - keep extending this locus
      $current_locus_stop = $stop if $stop > $current_locus_stop ;
    } # if gap > threshold
          
    # Add counts to current locus (can be a new one)
    $self->add_to_locus($$loci_ref{$locus_num}, $$hits_ref{$hit} , $samples_ref) ;
          
  } # foreach hit
  
  # The last locus after the last hit still
  # needs to be checked 
  # TODO  code duplication (this is the same as in the loop above) -> make sub
  if ($$loci_ref{$locus_num}{total_weighted_count} && $$loci_ref{$locus_num}{total_weighted_count} >= $min_hits ) { # keep it
    my $locus_pos = $current_chrom.';'.$current_locus_start .';'.$current_locus_stop ;
    $$loci_ref{$locus_num}{position} = $locus_pos ;
  } else { # not above threshold - delete 
    delete $$loci_ref{$locus_num} ;
  }	# if hit > threshold
} # process loci 
		
		
	
#####################################
# add_to_locus
# Called from generate_loci to add a hit
# to a locus. 
# Receives:
# hit_data->{ match_count => m,
#  read_count => { sample1 => r1,
#  sample2 => r2 } }	
# Normalised count:
# First weight the number of reads by the number
# of matches to the genome, i.e. 100 reads at this
# position which also match to 100 other positions
# give a count of 1.
# Then normalise to the number of genome matching reads
# in the sample multiplied by a factor (1000)
# The final score is the -repetition-weighted count
# per 1000 genome-matching reads for a given sample
sub add_to_locus{
  my $self = shift ;
  my $locus_ref = shift ;
  my $hit_data_ref = shift ;
  my $samples_ref = shift ;
  
  # each hit refers to one sequence
  # and therefore has a count of
  # matches of this seq to the genome
  # This only depends on the sequence
  # and is therefore sample-independent
  my $match_count = $$hit_data_ref{match_count} ;
  
  foreach my $sample (keys %{$$hit_data_ref{read_counts} } ) {
    my $read_count = $$hit_data_ref{read_counts}{$sample} ;
    my $weighted_count = $read_count / $match_count ;
    
    my $total_genome_matching_reads = $$samples_ref{$sample}{total_matching_reads} ;
    
  # Normlise: weigthed number of reads per 1 million genome matching reads
    my $normalised_count = $self->norm_count($weighted_count, $total_genome_matching_reads) ;
    
    $$locus_ref{read_counts}{$sample} +=$read_count ;
    $$locus_ref{weighted_counts}{$sample} +=$weighted_count ;
    $$locus_ref{normalised_weighted_counts}{$sample} += $normalised_count ;
    $$locus_ref{total_weighted_count} += $weighted_count ;
    $$locus_ref{num_unique}{$sample} += $read_count if $match_count == 1 ;  
  } # foreach sample

} # add_to_locus
	
  
#############################
# normalise counts:
# devide by total number of matching
# reads and multiply by a factor
# 1 million so that we return
# the normalised count in a unit:
# "count per 1 million matching reads" 
sub norm_count{
  my $self = shift;
  my $count = shift ;
  my $total_matching = shift ; # total number of genome matching reads
  return 0 unless $total_matching ;
  return ( ($count/$total_matching) * 1000000 );
 }

	
#############################
# sort_by_chrom_and_start
# input: match positions in 
#        format "chrom_id;start;stop;strand"
# Sorting order: by chromosome and then
# by start position.
# Chromosomes are sorted by number or
# ascii if they are strings
# This will work for both, hits and locus positions
# (strand is ignored)
sub sort_by_chrom_and_start{
  my ($achrom, $astart) = ($a=~/^(.+?);(\d+)/) ;
  my ($bchrom, $bstart) = ($b=~/^(.+?);(\d+)/) ;
  die "problem in sorting routine\n" unless $achrom && $bchrom && $astart && $bstart ;
  
  my $chrom_comp ;
  
  # We assume that the first character decides the format
  if ($achrom=~/^\d+/ && $bchrom=~/^\D+/ ) { # $achrom is a number and $bchrom a letter
    $chrom_comp =-1 ; # number always before letter
  } elsif ( $achrom=~/^\D+/ && $bchrom=~/^\d+/ ) { # other way around
    $chrom_comp =1 ; 
  } elsif ( $achrom=~/^\d+/ && $bchrom=~/^\d+/ ) { # both are numbers
    $chrom_comp = $achrom <=> $bchrom ; # compare numbers
  } elsif ( $achrom=~/^\D+/ && $bchrom=~/^\D+/ ) { # both are letters
    $chrom_comp = $achrom cmp $bchrom ; # compare ascii
  }
  
  # Now we have -1 or +1 for chromosome ID and we can return that
  # or a zero if they are the same -> then compare start positions
  return $chrom_comp || $astart <=> $bstart ;
} #



#############################
# compare_loci
# Get the comparison measures:
# log2 ratio (log2A - log2B)
# and average expression
# (also log2 average but
# this gnerates some unintuitive
# results on the normalised counts
# which can be <<1)
#############################
sub compare_loci{
  my $self = shift;
  my $loci_ref = shift;
  my $samples_ref = shift ;
  my $pseudocount = shift ;
  
  # Get an ordered list of sample names
  my @sample_names ; 
  foreach my $sample (sort keys %{$samples_ref}) {
    push @sample_names , $sample ;
  }
	
  # The pseudocount can be set as as a command line
  # parameter. It is zero by default
  # We must normalise this count too otherwise it
  # would skew the anaylsis if the two samples
  # were of different size
  my $norm_pseudocount_s1 =  $self->norm_count($pseudocount, $$samples_ref{$sample_names[0]}{total_matching_reads} ) ;
  my $norm_pseudocount_s2 =  $self->norm_count($pseudocount, $$samples_ref{$sample_names[1]}{total_matching_reads} ) ;
  
  # Traverse all loci - we are not printing so no need to sort here
  foreach my $locus_num (keys %{$loci_ref}) {
    # Calculations
    # Add pseudocounts to all counts (this can be zero)   
    my $norm_weighted_count_s1 =   $$loci_ref{$locus_num}{normalised_weighted_counts}{$sample_names[0]} || 0 ;
    my $norm_weighted_count_s2 =   $$loci_ref{$locus_num}{normalised_weighted_counts}{$sample_names[1]} || 0 ; 
    my $sample1_norm_w_count = $norm_weighted_count_s1 + $norm_pseudocount_s1;
    my $sample2_norm_w_count = $norm_weighted_count_s2 + $norm_pseudocount_s2; 
    $$loci_ref{$locus_num}{log2_ratio} = $self->log2_ratio($sample1_norm_w_count,$sample2_norm_w_count) ;
		
    # log2 avg expression: A = 0.5 * (log2 SAMPLE1 + log2 SAMPLE2)  
    # This can not be calculated if one of the counts is zero (return "n/a" in that case) 
    $$loci_ref{$locus_num}{log_avg_expr} = $self->log_avg_expr($sample1_norm_w_count, $sample2_norm_w_count) ;
    
    # Problem with log average: 
    # In some cases, the normalised counts can be very low. The avg log2 is then often negative
    # which is unintuitive.
    # This is the simple average expression level
    $$loci_ref{$locus_num}{avg_expr} = 0.5 * ($sample1_norm_w_count + $sample2_norm_w_count) ;
		
  } # each locus
  

} # compare_loci

    
#############################
# rank_loci
# Traverse all loci after all
# counts were done in
# process_loci and get two ranks
# for each: one for the log2 ratio
# and one for the average expression
# 
sub rank_loci{
  my $self = shift ;
  my $loci_ref = shift ;
  my $locus_ranks_ref = shift ;
  my $samples_ref = shift ;
  my $ratio_rank_weight = shift ;
  my $avg_rank_weight = shift ;
          
  ## If no pseudocount is used in combination with
  ## log2 averages we will have some "n/a" values
  ## which can not be ranked. Alternatively, we could
  # # then use the maximum expression rather than avg	
          
  # Traverse all loci
  # Sort by highest to lowest absolute log2 ratio
  # with "+/-inf" always first (highest regulation)
  # The "best" loci will have the lowest rank number.
  # Share rank if log2_ratio equal
  my $rank_number = 0 ;
  my $last_value = '' ;
  #foreach my $locus_num (sort log2_ratio_sort keys %{$loci_ref} ) {
  # see commented out subroutine for explanation of the sorting
  foreach my $locus_num (sort { 
    my $aratio = $$loci_ref{$a}{log2_ratio} || 0 ;                                                                                      my $bratio = $$loci_ref{$b}{log2_ratio} || 0 ;
    return 0 if $aratio=~/inf/ && $bratio=~/inf/ ;
    return -1 if $aratio=~/inf/ && $bratio!~/inf/ ;
    return 1 if $aratio!~/inf/ && $bratio=~/inf/ ;
    return abs($bratio) <=> abs($aratio) ; } 
    keys %{$loci_ref} ) {
    
    my $current_value = $$loci_ref{$locus_num}{log2_ratio} ;
          
    # If this is not the same log2_ratio as last time then
    # we can assign a new rank number to the next locus,
    # otherwise the next locus will share this rank
    unless ($current_value eq $last_value) {
            ++$rank_number ;
            $last_value = $current_value ;
    }
    
    $$loci_ref{$locus_num}{ratio_rank} = $rank_number ;
          
            
  } # each locus

  # Now sort by expression level and assign ranks for that too
  # The lower the rank number the better (higher avg expression)
  # Use "simple" averages for now, not log2 averages
  # TODO: decide which avg to use
  # At this stage we can also assign the weighted rank sum score
  # calculated as follows:
  # rs = rw * rr + aw * ar
  # rs = rank sum, rw = ratio weight, rr = ratio rank,
  # aw = averages weight, ar = averages rank
  # so those that have both high regulation and high
  # average expression get a low sum-rank and are the
  # best candidates
  # The weights for the sum rank are 1 by default
  $rank_number = 0 ;
  $last_value = '' ;
  #foreach my $locus_num (sort avg_expr_sort keys %{$loci_ref} ) {
  # see commented out subroutine for explanation of the sorting
  foreach my $locus_num (sort { 
    my $aavg = $$loci_ref{$a}{avg_expr} || 'n/a' ;
    my $bavg = $$loci_ref{$b}{avg_expr} || 'n/a' ;
    return 1 if $aavg!~/^-?\d/ && $bavg=~/^-?\d/ ;
    return -1 if $bavg!~/^-?\d/ && $aavg=~/^-?\d/ ;
    return 0 if $aavg!~/^-?\d/ && $bavg!~/^-?\d/ ;
    return $bavg <=> $aavg ;
    } keys %{$loci_ref} ) {
          
    my $current_value = $$loci_ref{$locus_num}{avg_expr} ;
          
    # If this is not the same log2_ratio as last time then
    # we can assign a new rank number to the next locus,
    # otherwise the next locus will share this rank
    unless ($current_value eq $last_value) {
            ++$rank_number ;
            $last_value = $current_value ;
    }
    
    $$loci_ref{$locus_num}{avg_rank} = $rank_number ;
    
    # We can now calculate the (weighted) rank sum
    # By default, both ranks (ratio and average) are weighted
    # equaly 
    $$loci_ref{$locus_num}{rank_sum} = $ratio_rank_weight * $$loci_ref{$locus_num}{ratio_rank} + $avg_rank_weight * $$loci_ref{$locus_num}{avg_rank} ; 
  } # each 



} # rank_loci
	
#############################
# generate_result_files
# in working_dir/results
sub generate_result_files{
  my $self = shift ;
  my $loci_ref = shift ;
  my $working_dir = shift ;
  my $samples_ref = shift ; # sample short names (usually S1,S2)
  my $output_options_ref = shift ; 
  
  # We extract this from the job cache where
  # we have a mapping of the upload file name to
  # sample name and a mapping of short sample names
  # (S1,S2) to the file names (this was done in two steps
  # to prevent the property filter module from having to 
  # handle sample names, which are only used for this tool)
  my %sample_short_to_full_name ;
  
  # Before starting to write, make sure working_dir
  # is still there (last check could have been
  # more than an hour ago)
  die " could not find/read working directory\n" unless -d $working_dir ;

  # Name of the locus file:
  # loci_SAMPLE1_SAMPLE2.csv
  my $result_file_name = "loci" ;
  foreach my $sample (sort keys %{$samples_ref} ) {
    my $sample_file = $self->job->cache->{sample_to_filenames}{$sample};
    my $sample_full_name = $self->job->cache->{file_to_sample_param}{$sample_file};
    $sample_short_to_full_name{$sample} = $sample_full_name ;
    $result_file_name .='_'.$sample_full_name;
  }
  $result_file_name .='.csv' ;
	
  # Open result files and print header
  my $resultfile = $working_dir . '/results/'.$result_file_name ;
  open (RESULT_FILE, ">", $resultfile) || die " could not open result file\n" ; 
	
  # Print some information into the first rows
  print RESULT_FILE "\"sRNA locus comparison" ;
  print RESULT_FILE " (single sample mode)" if $$output_options_ref{use_single_sample} ;
  print RESULT_FILE "\"\n" ;
  
  print RESULT_FILE "\"Weighted count: Each read count is divided by the number of its matches to the genome\"\n" ;
  print RESULT_FILE "\"Normalised count: Sum of weighted counts normalised to total number of reads matching the genome in each dataset. Given in counts per 1 million genome-matching reads\"\n" ;
  print RESULT_FILE "\"PLEASE NOTE: if a pseudocount was used, the log-ratio and averages will be calculated using normalised count + pseudocount\"\n" unless $$output_options_ref{use_single_sample};
  # Get the display name of the genome used
  print RESULT_FILE "\"Genome used: ".$self->job->get_genome_dname."\"\n" ;
  
  # Add legend for the samples and the stats
  # The number of reads/valid reads, nr sequences
  # is opbtained from the config file. The number 
  # of genome-matching reads is obtained in this script
  # Put sample names into an array so we freeze the order 
  # in which we print their data
  print RESULT_FILE "\n\"Samples:\"\n" ;
  my @ordered_sample_names ; 
  my $i = 1 ;
  foreach my $sample (sort keys %{$samples_ref}) {
    # sample names used here are "shorthand" (S1,S2),
    # get full name for this sample from config
    my $sample_full_name = $sample_short_to_full_name{$sample} ;
    die "Could not read full sample name from job cache. Please contact administrator\n" unless $sample_full_name ;
    
    print RESULT_FILE "\"$sample = $sample_full_name\"\n" ;
    push @ordered_sample_names , $sample ;
    ++$i;
  }
  
  # Print read count table
  print RESULT_FILE "\"Read counts\"\n";
  print RESULT_FILE "\"\"" ;
  foreach my $sample (sort keys %{$samples_ref}) {
    print RESULT_FILE ",\"$sample total\",\"$sample non-redundant\"";
  }
  print RESULT_FILE "\n" ;
  
  # stages are "input", "filtered: sequence properties", "filtered: t/rRNA"
  foreach my $stage ( @{$self->job->cache->{read_counts}} ) {
    my ($stage_name) = (keys %$stage) ; # there is only one key 
    print RESULT_FILE "\"$stage_name\"" ;
    foreach my $sample (@ordered_sample_names) { 
      print RESULT_FILE ','.$$stage{$stage_name}{$sample}{total};
      print RESULT_FILE ','.$$stage{$stage_name}{$sample}{nr} ;
    }
    print RESULT_FILE "\n" ;
	}
  # Add genome-matching counts
  # generated in this script
  print RESULT_FILE "\"genome-matching\"" ;
  foreach my $sample (@ordered_sample_names) { 
    print RESULT_FILE ','.$$samples_ref{$sample}{total_matching_reads} ;
    print RESULT_FILE ','. $$samples_ref{$sample}{nr_matching_reads};
  }
  print RESULT_FILE "\n\n" ;
  
  # Header for the data
  # Locus position
  print RESULT_FILE "\"locus Number\",\"chromosome\",\"start position\",\"end position\",\"length\"" ;
  # Counts
  foreach ( @ordered_sample_names) {
    print RESULT_FILE ",\"raw count $_\",\"Weighted count $_\",\"Normalised count $_ \"";
    print RESULT_FILE ",\"Uniquely matching reads $_\"" if $$output_options_ref{show_unique} ;
	}
  # log ratios and average expression
  # (not in single sample mode)
  unless ($$output_options_ref{use_single_sample}) {
    print RESULT_FILE ",\"log2 ratio\"";
    print RESULT_FILE ",\"avg norm. count\"";
    print RESULT_FILE ",\"log2-ratio rank\"";
    print RESULT_FILE ",\"avg-based rank\"";
    print RESULT_FILE ",\"weighted rank sum \"";
  }
	# ..and the hyperlinks to genome browsers if any
  foreach my $gb_link ( sort @{$$output_options_ref{gb_links}}) {
    print RESULT_FILE ",\"link to $gb_link genome browser\"" ;
  }
  print RESULT_FILE "\n" ; # end header
	
  
  # Traverse loci
  # The locus numbers were assigned in chromosome / start-pos
  # order, so traversing by locus number also sorts by position
  # and print csv format
  foreach my $locus_num (sort {$a <=> $b} keys %{$loci_ref}) {
  
    my $pos = $$loci_ref{$locus_num}{position} ;
    my ($chrom,$start, $stop) ;
    if ($pos=~/^(.+?);(\d+);(\d+)$/ ) {
      ($chrom,$start, $stop) = ($1, $2, $3) ;
    } else {
      die " could not parse position string $pos\n" ;
    }
    my $length = $stop - $start + 1 ;
    
    # Format the counts: use scientific notation for very high or low numbers
    my %sample_read_count ;
    my %sample_w_count ;
    my %sample_norm_w_count ;
    my %sample_unique_count ;
    my $log_ratio ;
    my $log_avg_expr ;
    my $avg_expr ;
    my $ratio_rank ;
    my $avg_rank  ;
    my $rank_sum ;
    foreach my $sample_name (@ordered_sample_names) { 
      $sample_read_count{$sample_name} = $self->format_counts( $$loci_ref{$locus_num}{read_counts}{$sample_name} );
      $sample_w_count{$sample_name} = $self->format_counts( $$loci_ref{$locus_num}{weighted_counts}{$sample_name} );
      $sample_norm_w_count{$sample_name} = $self->format_counts( $$loci_ref{$locus_num}{normalised_weighted_counts}{$sample_name} ) ;
      $sample_unique_count{$sample_name} = $self->format_counts( $$loci_ref{$locus_num}{num_unique}{$sample_name} );
    }
          
    unless ($$output_options_ref{use_single_sample}) {
      $log_ratio = $self->format_counts( $$loci_ref{$locus_num}{log2_ratio} );
      $log_avg_expr = $self->format_counts( $$loci_ref{$locus_num}{log_avg_expr} ) ;
      $avg_expr = $self->format_counts( $$loci_ref{$locus_num}{avg_expr} );
      $ratio_rank = $$loci_ref{$locus_num}{ratio_rank} ;
      $avg_rank = $$loci_ref{$locus_num}{avg_rank} ;
      $rank_sum =  $$loci_ref{$locus_num}{rank_sum} ;
    }

    # Print results in csv format
    # log ratio and avg expression can be text in some cases, so must be quoted
    print RESULT_FILE "$locus_num,\"$chrom\",$start,$stop,$length" ;
    foreach my $sample_name (@ordered_sample_names) { 
      print RESULT_FILE ",$sample_read_count{$sample_name},$sample_w_count{$sample_name},$sample_norm_w_count{$sample_name}" ;
      print RESULT_FILE ",$sample_unique_count{$sample_name}" if $$output_options_ref{show_unique} ;
    }
    unless ($$output_options_ref{use_single_sample}) {
      print RESULT_FILE ",\"$log_ratio\",\"$avg_expr\",$ratio_rank,$avg_rank,$rank_sum";
    }
    # Add links to genome browser if requested.   
    # In the config, we have templates for URLs and URLS parameters
    # for the genomebrowsers that we put together here
    foreach my $gb_link ( sort @{$$output_options_ref{gb_links}}) {
      my $base_url = $self->job->cfg->{gblinks}{$gb_link}{base_url};
      my $params = $self->job->cfg->{gblinks}{$gb_link}{params} ;
      my $link = '--';
      if (defined $base_url || defined $params) {
        $link = $self->generate_gb_link($chrom,$start,$stop, $base_url, $params  ) ;
      }
      print RESULT_FILE ",\"$link\"" ;
    }
          
    print RESULT_FILE "\n" ;

  } # each pos

  close RESULT_FILE ;
} # generate result files
	


######################
# log average expression
# A = 0.5 * (log2 SAMPLE1 + log2 SAMPLE2)
# as used in microarray analysis
# This can not be calcualted if
# one of the samples is zero and no
# pseudocount is added
# -> return "n/a"
sub log_avg_expr{
  my $self = shift ;
  my $sample1 = shift ;
  my $sample2 = shift ;
  
  return 'n/a' unless $sample1 && $sample2 ;
  
  return ( 0.5 * ($self->log_n($sample1,2) + $self->log_n($sample2,2)) );
}


#######################
# Calculate log2 ratio
# of two expression levels
# (normalised weighted counts)
# as log2(A/B) (=log2A - log2B)
# return +inf if A=0 and
# -inf if B=0.
# Perl log returns the
# log of base e. Devide
# by loge(2) to get
# log2(n)
sub log2_ratio{
  my $self = shift ;
  my $a = shift ;
  my $b = shift ;
  
  # A locus always has a non-zero count
  # for one of the samples
  if (!$a) {
    return '+inf';
  } elsif (!$b ) {
    return '-inf';
  } else {
    # get log2 difference for reads
    return ($self->log_n($a, 2) - $self->log_n($b, 2) );
  }
} # log_ratio

########################
# Perl log returns the
# log of base e. Devide
# by loge(n) to get
# logn(n), e.g. for n=2
sub log_n{
  my $self = shift ;
  my $n = shift ;
  my $base = shift ;
  return (log($n)/log($base));
}
	
# Create link to tair grbowse from a position
# The params line contains placeholders in the
# format *START*, *STOP*, *CHROM* which we 
# replace with the required values here
sub generate_gb_link {
  my $self = shift;
  my $chrom = shift ;
  my $start = shift ;
  my $stop = shift ;
  my $base_url = shift ; # genome browser URL w/o the parameters
  my $params = shift ; # parameter part of the URL (?....)
  
  # asrp doesn't have mito and chloroplast
  return 'n/a' if ($base_url =~/asrp/ && $chrom=~/^[MC]$/i);
  
  $params=~s/\*START\*/$start/ ;
  $params=~s/\*STOP\*/$stop/ ;
  $params=~s/\*CHROM\*/$chrom/ ;
	
  return $self->excel_link($base_url . $params, undef) ;
}

######################
# Format numbers for 
# tables.
# Switch to scientific
# notation if above or below
# thresholds
sub format_counts{
  my $self = shift ;
  my $n = shift ;
  return 0 unless $n ; # leave zeros and turn undef into zero
  return $n unless $n=~/^-?\d+/; # only apply to numbers
  if (abs($n)>9999 || abs($n)<0.01) {
    return sprintf ("%.2e", $n) ;
  } else {
    return sprintf ("%.2f", $n) ;
  }
	
} # format_counts


# Generate a link for use in Excel
# suitbale for csv output, i.e. with 
# duplicated quotes
sub excel_link {
  my $self = shift;
  my $href = shift ;
  my $text = shift ;
	
  # Not so nice but works in Open Office too
  # if we omit the link text
  # Excel and OOcalc differ in the way
  # they separate ref and display text
  # (, for Excel ; for OO)  
  return "=HYPERLINK(\"\"$href\"\")" ;	
}


1;