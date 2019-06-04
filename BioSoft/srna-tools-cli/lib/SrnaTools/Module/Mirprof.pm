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
# It builds a table of expression levels for known miRNAs.
# Matches are read from Patman output supplied to this script.
# The script first builds a hash of sequences with arrays of their
# matching miRNAs, then it transposes these to a hash keyed by
# the match strings, which can be combined by deleting organism
# suffixes, family member information and number of mismatches.
# 
#
# Arguments
# --working_dir             Working-directory for this project
# --proj                    Name of job for output
# -p FILE                   Name of Patman result file (expected in
#                           WORKING_DIR/data/)
# --ignore_org              ignore organism in miRNA ID for grouping
# --ignore_mis              ignore mismatches in match string for grouping
# --collapse_match_groups   combine matches by concatenating all match IDs
#                           of each sequence and using the concatenated
#                           strings as IDs. Each of the resulting groups
#                           only has unique matches (total and weighted
#                           match counts are always equal) 
# --group_family            Group family members together (ignore a,b,c suffix)
# --group_variant           Group different variants (everything after the family
#                           suffix in miRNA ID, e.g. -1, .1, -5p) together
# --keep_best               Keep only "best matches" for each sequence; e.g.
#                           if a sequence has matches with zero mismatches then
#                           all matches with 1 or more mismatches are discarded
# --alt_bin                 Optional alternative path to the zip binary
#                           This will normally be lib/bin and the executable
#                           should be copied or symlinked in that directory.

#
# We generate the output in the subdir "results" in the working directory. 

package SrnaTools::Module::Mirprof ;
use base SrnaTools::Module ;
use strict ;
use warnings;

sub run{
my $self = shift ;

#############################
# Parameters and
# declarations
#############################
my $module_name = "miRNA profiling" ; # for error log
my $infile = $self->param('infile') ;
my $infile_full_path ;
my $ignore_org = $self->param('ignore_org') ;
my $ignore_mis = $self->param('ignore_mis');
my $collapse_match_groups = $self->param('collapse_match_groups');
my $group_family = $self->param('group_family');
my $group_variant = $self->param('group_variant');
my $keep_best = $self->param('keep_best');
my %output_opts;
my %srna_matches ;
my %srnas ;
my $job_name = $self->job->job_name;
my $working_dir = $self->job->job_working_dir ;
my $data_dir ;

#############################
# Check parameters and files
#############################

eval {
  die "missing working directory parameter - please contact administrator\n" unless $working_dir ;
  die "missing job name parameter - please contact administrator\n" unless $job_name ;
  
  $data_dir = $working_dir.'/data/' ;
  $infile_full_path = $working_dir.'/data/'.$infile ;
  
  # Check files and directories
  die "could not find/read working directory - please contact administrator\n" unless -d $working_dir ;
  die "data directory not found in working directory - please contact administrator\n" unless -d $data_dir ;  
  die "could not find/read Patman result file - please contact administrator\n" unless -r $infile_full_path ;
  
  die "could not find zip executable - please contact administrator\n" unless $self->binary_in_path('zip') ;
  
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error in parameters: $@") ;
}

# Set some output options
# TODO get rid of this
$output_opts{ignore_org} = $ignore_org ;
$output_opts{ignore_mis} = $ignore_mis ;
$output_opts{collapse_match_groups} = $collapse_match_groups ;
$output_opts{group_family} = $group_family ;
$output_opts{group_variant} = $group_variant ;
$output_opts{keep_best} = $keep_best ;
$output_opts{job_name} = $job_name ;

if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error reading config file") ;
}

# update status on server
$self->update_status("$module_name: generating miRProf output") ;

#############################
# Parse hits from Patman.
# Each sequence is stored with
# a read count and a hash of
# match-strings in the format
# miRNA_ID(mismatches).
#############################

eval{ # Catch errors for this part
  $self->parse_hits( 
    \%srnas,
    $infile_full_path,
    \%output_opts, # for ignore_org, ignore_mis
  ) ; 
} ;
# Do this if we had an exception
if($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while reading Patman result file: $@") ;
}

#############################
# The srnas hash collects match 
# (miRNA) IDs by sequence - we 
# now need to transpose this into
# a hash keyed by the concatenated
# match IDs. 
#############################
eval{ # Catch errors for this part
  $self->transpose_result_hash( 
    \%srnas,
    \%srna_matches,
    \%output_opts,
  ) ;
} ;
# Do this if we had an exception
if($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while processing data: $@") ;
}

#############################
# Generate output
#############################
eval{ # Catch errors for this part
  my $dest_dir = $working_dir.'/results/mirprof_results' ;
  mkdir $dest_dir or die "could not generate result directory\n" ;
  $self->generate_table(
    \%srna_matches,
    \%output_opts,
    $dest_dir,
  ) ;
  $self->generate_fasta(
    \%srna_matches,
    \%output_opts,
    $dest_dir,
  ) ;
  # zip results: -q quite, -r recursive, -j exclude full path to dir
  my $cmd = "zip -qrj $dest_dir $dest_dir" ;
  system($cmd) == 0 or die "Failed to compress result files\n" ; 
  system("rm -rf $dest_dir") ; # remove original directory
} ;
# Do this if we had an exception
if($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while generating result files: $@") ;
}
return 1;
} #run


#############################
# parse_hit
# Parse the output from patman
# where abundances of each match
# sequence were stored in the
# ID line in the format
# SEQUENCE_sample1name(abundance)[_sample2name(abundance)
#############################
sub parse_hits{
  my $self = shift ;
  my $srnas_ref = shift ;
  my $pat_out_file = shift ;
  my $output_opts_ref = shift ; 
  
  # Open the patman output file
  open(PAT, $pat_out_file) or die "could not open PatMaN output file\n" ;
  while ( my $patline = <PAT>) {
    
    # parse the hit
    # patman stores the entire ID line from
    # FASTA, not just ">\S+" like Bioperl
    my ($match_id, $seq, $seq_abund, $mismatches) ;
    if ($patline=~/^(\S+).*?\t([AGCT]+)(_\S+)\t\d+\t\d+\t[+-]\t(\d+)/  ) {
      ($match_id, $seq, $seq_abund, $mismatches) = ($1,$2,$3, $4) ;
    } else {
      die "Patman result line $patline could not be parsed - please contact administrator\n" ;
    }
    
    # Add match (miRNA) ID to list of matches 
    # for this sequence. Add number of mismatches
    # in the format miRNAid(MISMATCHES)
    my $match_string = $match_id.'('.$mismatches.')' ;
    my ($seq_count) = ($seq_abund=~/_\S+?\((\d+)\)/) ;
    $$srnas_ref{$seq}{count} = $seq_count;
    push @{$$srnas_ref{$seq}{matches}}, $match_string ;
  } # parser (patman file)
 close PAT ;     
} # parse_hits


#############################
# transpose result hash:
# build a new hash keyed
# by matches rather than sRNAs.
# This is required to print the
# results by match (or match group)
# in table and FASTA file.
# Grouping and filtering for best
# hits only is done at this step
# (it would not be possible to filter
# for best matches before all sRNAs and
# their matches have been read from
# the patman file, which is why the
# 2 datastructures are needed).
#############################
sub transpose_result_hash{
  my $self = shift ;
  my $srnas_ref = shift ;
  my $srna_matches_ref = shift ;
  my $output_opts_ref= shift ;
  
  # traverse sRNAs 
  foreach my $seq (keys %$srnas_ref) {
    my $seq_count = $$srnas_ref{$seq}{count} ; # Get abundance of sRNA
    
    # If the keep_best option is chosen,
    # discard all miRNAs that have more
    # mismatches than the best hit(s)
    $self->discard_worse_matches($$srnas_ref{$seq}{matches}) if $$output_opts_ref{keep_best};
    
    # Group matches according to options. 
    my %grouped_matches = $self->do_group_matches($$srnas_ref{$seq}{matches}, $output_opts_ref) ;
    
    # Collapse into org="all combined" and only one 
    # match string if collapse_match_groups (match signature) used.
    # The match string is simply all matches concatamerised into one
    # semicolon delimited string, e.g. an sRNA mathing mir319a and mir319b
    # would have the signature mir319a; mir319b, to which it uniquely matches
    $self->do_collapse_match_groups(\%grouped_matches) if $$output_opts_ref{collapse_match_groups};
    
    # Process matches separatedly by organisms
    # (if ignore_org used this is only one: "all combined")
    foreach my $org (sort keys %grouped_matches) {
      my $num_matches = keys %{$grouped_matches{$org}} ; # number of matching miRNAs in this organism
    
      # Now add all (combined) matches for this sequence
      # to the global list of matches with counts
      foreach my $match (sort keys %{$grouped_matches{$org}} ) { 
        $$srna_matches_ref{$org}{$match}{raw_count} += $seq_count ;
        $$srna_matches_ref{$org}{$match}{w_count} += $seq_count / $num_matches ;
        # Keep the sequence for fasta output (with its count)
        $$srna_matches_ref{$org}{$match}{seqs}{$seq} = $seq_count ;
      }
      # don't need the hash element keyed by sequence anymore
      delete $$srnas_ref{$seq} ;
    } # each org
  } # each seq
} #transpose_result_hash


#############################
# Combine matches according
# to options. Keep miRNAs from
# different organisms separate because
# we don't want to split up counts across 
# organisms (e.g. a sequence matching ath-miR156
# and osa-miR156 shouldn't get a weigthed 
# count of 0.5). If organisms are ignored, 
# we add all matches to organism 'all combined'
#############################
sub do_group_matches{
  my $self = shift ;
  my $matches_ref = shift ;
  my $output_opts_ref = shift ;
  
  my %grouped_matches ;

  foreach my $match_string ( @$matches_ref ) {
        
    # capture organism, family and variant information
    # from the miRNA ID
    $match_string=~/
      ^(\w+?)-             # organism acronym, ath, cre etc. - we need this later
      (miR-?)?             # in almost all cases the name starts with "miR-" (animals)
                           # or "miR" (plants) but dme-bantam or cel-let-7 
                           # are also valid as is miR-iab. 
      (\d+|\w+(?:-\d+)?)   # main miRNA ID
                           # plants: miR123
                           # animals: miR-123 or let,bantam or miR-iab
                           # astring name can be followed by a number (let-7)
                           # but don't capture that number separatedly
      ([a-z])?             # Family member information (if any), examples:
                           # plants: miR162a
                           # animals: let-7a, miR-29a
      ([^(]*)              # Anything between the ID (family member)
                           # and the number of mismatches, starting with "("
                           # is the "variant", 
                           # e.g. miR-125a-5p is the 5' mature sequence
                           # mir161.1 and 161.2 are overlapping mature seq 
                           # from same hairpin
                           # and miR219-1-3p and miR219-2-3p are the same mature 
                           # seqs from two distinct precursors  
      (\(\d+\))?           # number of mismatches to sRNA in brackets (from miRProf)
    /ix ;
      
    # some parts may be empty strings
    my ($org, $mir_str, $id, $fam_mem, $variant, $mismat) = ($1, $2||'', $3, $4||'', $5||'', $6||'') ;
    
    # reconstruct the final match string using
    # only the parts as requested by options
    my $final_match_str = $mir_str . $id ;
    $final_match_str .= $fam_mem unless $$output_opts_ref{group_family} ;
    $final_match_str .= $variant unless $$output_opts_ref{group_variant} ;
    $final_match_str .= $mismat unless $$output_opts_ref{ignore_mis} ;
    
    # Either store by organism or combine all into one
    $org = $$output_opts_ref{ignore_org} ? 'all combined' : $org ;
    $grouped_matches{$org}{$final_match_str} = 1 ;

  }
  return %grouped_matches ;
} # do_group_matches


#############################
# If we are combining by "match signature", we only have one key
# per sequence (the concatanation of all miRNA IDs after all grouping
# is done) per organism 
#############################
sub do_collapse_match_groups{
  my $self = shift ;
  my $grouped_matches_ref = shift ;
  
  foreach my $org (keys %$grouped_matches_ref) {
    my $match_string = join('; ', sort keys %{$$grouped_matches_ref{$org}}) ;
    # replace original match strings with this new single
    # concatanated string 
    $$grouped_matches_ref{$org} = {$match_string => 1} ; 
  }
  
} # get_collapsed_match_groups

#############################
# discard match-IDs that have
# more mismatches than the 
# best one(s)
# At the moment, this will
# discard matches from all
# species. If matches of each
# sequence would be stored by
# organism, we could change this
# to discard only within same organism
#############################
sub discard_worse_matches{
  my $self = shift ;
  my $match_ref = shift ;
	
  return if @$match_ref < 2 ;
  
  my @temp ;
  my $best;
  my $is_first = 1 ;
  
  # Match strings contain the number of
  # mismatches in brackets. This is done
  # before the grouping options are applied,
  # so the mismatch number must be present)
  foreach my $match_string (sort sort_by_mismatch @$match_ref ) {
    my ($num_mismatch) = ($match_string=~/\((\d+)\)/) ;
    die "Could not extract number of mismatches" unless defined $num_mismatch ;
    if ($is_first) {
      $best = $num_mismatch ;
      $is_first = 0 ;
    }
    last if $num_mismatch > $best ;
    push @temp, $match_string ;
  }
  # re-assign filtered list to original array
  @$match_ref = @temp ;
  return ;
} #discard_worse_matches


#############################
# sort by mismatch, which is 
# stored in match string as
# ID(MISMATCH)
#############################
sub sort_by_mismatch{
  my ($aa) = ($a=~/\((\d+)\)/) ;
  my ($bb) = ($b=~/\((\d+)\)/) ;
  return $aa <=> $bb ;
}


#############################
# Generate csv file
#############################
sub generate_table{
  my $self = shift;
  my $srna_matches_ref = shift ;
  my $output_opts_ref = shift ;
  my $dest_dir = shift ;
   
  # Get total read count from config
  my $total_reads = $self->get_total_reads_from_cache ;
  die "Could not get total read count\n" unless $total_reads;
  
  my $file_name = $$output_opts_ref{job_name}.'_profile.csv' ;
  open (RESULT , '>', $dest_dir.'/'.$file_name) or die "Could not generate result file\n" ;
  
  print RESULT "\"miRProf results for $$output_opts_ref{job_name}\"\n" ;
  print RESULT "\"Normalised count: count of matching sequence reads normalised to total number of reads after last filtering step (see table below). Given in parts per million\"\n" ;
  if (defined $self->job->param('genome') ) {
    print RESULT "\"Genome used for filtering: ".$self->job->get_genome_dname."\"\n" ;
  }
  
  # print name of mirbase database used
  print RESULT "\"miRBase database: ".$self->job->param('mirbase_db') ;
  print RESULT ", number of allowed mismatches: ".$self->job->param('mismatches')."\"\n" ;
  # print used grouping options
  my $group_opts_string = $self->get_group_opts_string($output_opts_ref) ;
  print RESULT "\"Grouping options used: $group_opts_string\"\n" ;
  my $keep_string = $$output_opts_ref{keep_best} ? 'keep best matches only' : 'keep all matches';
  print RESULT "\"Mismatch filtering: $keep_string\"\n" ;
  print RESULT "\n" ;
  
  # print table of read counts
  # stages are "input", "filtered: sequence properties", "filtered: t/rRNA"
  # We only have one sample S1 at the moment
  print RESULT "\"Read counts\"\n";
  print RESULT "\"\",\"total\",\"non-redundant\"\n" ;
  foreach my $stage ( @{$self->job->cache->{read_counts}} ) {
    my ($stage_name) = (keys %$stage) ; # there is only one key 
    print RESULT "\"$stage_name\"" ;
    print RESULT ','.$$stage{$stage_name}{'S1'}{total};
    print RESULT ','.$$stage{$stage_name}{'S1'}{nr} ;
    print RESULT "\n" ;
  }
  
  # data
  # one block per organism (which could be
  # "all combined")
  foreach my $org ( sort keys %$srna_matches_ref ) {
    print RESULT "\n\"Organism: $org\"\n" ;
    # Header for the data
    print RESULT "\"matching miRNAs\",\"raw count\"" ;
    print RESULT ",\"weighted count\"" ;
    print RESULT ",\"normalised weighted count\"\n" ;
    foreach my $match ( sort keys %{$$srna_matches_ref{$org}} ) {
      print RESULT "\"$match\"" ;
      print RESULT ','.$$srna_matches_ref{$org}{$match}{raw_count} ;
      print RESULT ','.$self->format_counts($$srna_matches_ref{$org}{$match}{w_count}) ;
      
      # We normalise to total reads after last filtering step
      my $norm_count = ($$srna_matches_ref{$org}{$match}{w_count} / $total_reads) * 1e6 ;
      print RESULT ','.$self->format_counts($norm_count) ;
      print RESULT "\n" ;
    } # each match
  } # each org

  close RESULT ;
} # generate_table

######################
# format grouping options
# for table
######################
sub get_group_opts_string{
  my $self = shift ;
  my $output_opts_ref = shift ;
  my @opts_strings ;
  push (@opts_strings, 'match signature') if $$output_opts_ref{collapse_match_groups};
  push (@opts_strings, 'ignore mismatches') if $$output_opts_ref{ignore_mis} ;
  push (@opts_strings, 'combine organisms') if  $$output_opts_ref{ignore_org} ;
  push (@opts_strings, 'combine family members')  if $$output_opts_ref{group_family} ;
  push (@opts_strings, 'combine variants')  if $$output_opts_ref{group_variant} ;
  my $opts_string = @opts_strings ? join(', ', @opts_strings) : 'none' ;
  return $opts_string ;
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


#############################
# extract the total read
# count after the last filtering
# step (t/rRNA, genome)
#############################
sub get_total_reads_from_cache{
  my $self = shift ;
  my $last_filter = @{$self->job->cache->{read_counts}}[-1] ;
  my ($filter_name) = (keys %$last_filter) ; # there is only one key
  return $$last_filter{$filter_name}{'S1'}{total} ; # only one sample
  
} # get_total_reads_from_config


#############################
# Generate fasta file
#############################
sub generate_fasta{
  my $self = shift ;
  my $srna_matches_ref = shift ;
  my $output_opts_ref = shift ;
  my $dest_dir = shift ;
  
  my $file_name = $$output_opts_ref{job_name}.'_sequences.fasta' ;
  open (RESULT , '>', $dest_dir.'/'.$file_name) or die "Could not generate result file\n" ;
  
  foreach my $org ( sort keys %$srna_matches_ref ) {
    foreach my $match ( sort keys %{$$srna_matches_ref{$org}} ) {
      my $i = 0;
      # sort descending by raw count
      foreach my $seq (sort {$$srna_matches_ref{$org}{$match}{seqs}{$b} <=> $$srna_matches_ref{$org}{$match}{seqs}{$a} } keys %{$$srna_matches_ref{$org}{$match}{seqs}} ) {
        ++$i ;
        my $count = $$srna_matches_ref{$org}{$match}{seqs}{$seq} ;
        my $id = $org.'-'.$match.'_'.$i.'_'.$count.'x' ;
        $id=~s/\s/_/g ; # make sure we don't have spaces in ID (e.g. from "all combined")
        print RESULT ">$id\n$seq\n" ;
      } # each seq
    } # each match
  } # each org
  close RESULT ;
} # generate_fasta

1 ;