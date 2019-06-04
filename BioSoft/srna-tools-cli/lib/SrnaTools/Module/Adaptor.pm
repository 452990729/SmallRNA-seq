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
# It removes adaptor sequences from sRNAs provided in FASTA or
# FASTQ format and also filters by seq length
#
#####################################################################
package SrnaTools::Module::Adaptor ;
use base SrnaTools::Module ;
use strict ;
use warnings;

sub run{
my $self = shift ;

require Bio::SeqIO;

my $module_name = "Adaptor remover";
my $working_dir = $self->job->job_working_dir ;
my $infile = $working_dir.'/data/'.$self->param('infile');
my $minsize = $self->param('minsize');
my $maxsize = $self->param('maxsize');
my $adaptor_5 = $self->param('adaptor_sequence_5')||'';
$adaptor_5 =~ tr/a-z/A-Z/; # Make uppercase
$adaptor_5 =~ s/U/T/g; # Convert U's to T's

my $adaptor_3 = $self->param('adaptor_sequence_3');
$adaptor_3 =~ tr/a-z/A-Z/; # Make uppercase
$adaptor_3 =~ s/U/T/g; # Convert U's to T's

my $revcomp = $self->param('allow_revcomp'); # Do we search for reverse complement matches (0 = no, default. 1 = yes)

my $output_dir = $working_dir.'/results/'; # Where to write the file on the cluster

eval {
  die "missing working directory parameter - please contact administrator\n" unless $working_dir ;
  
  # Check files and directories
  die "could not find/read working directory - please contact administrator\n" unless -d $working_dir ;
  die "could not find/read input file - please contact administrator\n" unless -r $infile ;
  
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error in parameters: $@") ;
}


# Initialise counters:

my $total_count = 0; # total number of sequences in input file
my $processed_count_5 = 0; # number of sequences matching 5' adaptor
my $processed_count_3 = 0; # number of sequences matching 3' adaptor
my $processed_count_both = 0; # number of sequences matching both adaptors
my $processed_count_valid = 0; # number of sequences matching both adaptors and processed sRNA falls within the user specified size range

# Hash containing sequence length counts for stats:

my %lengths;

# update status on server
$self->update_status("$module_name: converting input sequences") ;

# The next step replaced the fastq2fasta.pl script
my $converted_input ;
eval {
  $converted_input = $self->convert_fastq_to_fasta($infile) ;
} ;
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while converting to FASTA format: $@") ;
}

$self->update_status("$module_name: extracting sRNAs") ;


eval {
my $seqs = Bio::SeqIO->new(-format => 'fasta', -file => "$converted_input") or SrnaTools::Exception::ModuleExecution->throw("$module_name: could not generate Bio::SeqIO object. Probably due to file not existing") ; # Read in sequences unsing BioPerl
open (FASTA_FINAL, ">", $output_dir.'/srnas_adapters_removed.fa') or die "could not open file for final FASTA output";

while (my $seq = $seqs->next_seq ) { # Get next sequence object
  $total_count++; # Increment counter for total number of sequences
  my $sequence = $seq->seq(); # Get the sequence string from the Seq object
  my $id = $seq->id(); # Get the sequence ID from the Seq object
  my $processed = 0; # May or may not have a 5' adaptor so only allow 3' matches if we match the 5' adaptor OR if 5' adaptor not specified
  if ($revcomp){ # If we do want to match the reverse complement then:
    my $rc = $seq->revcom->seq(); # Get the reverse complemented sequence string
    if (!$adaptor_5){ # If no 5' adaptor specified then:
      $processed++; # Increment this
    }
    elsif ($rc =~m/$adaptor_5(\w+)$/i){ # If there is a 5' adaptor specified and we match it then:
      $rc = $1; # capture the sequence after the adaptor match
      $processed_count_5++; # We've matched a 5' adaptor so increment this counter
      $processed++; # Ok, so we've matched the 5' adaptor so we'll keep track of this
      }
      if ($rc =~m/^(\w+)$adaptor_3\w*$/i){ # If we match the 3' adaptor then:
        my $processed_seq = $1; # Capture the sequence before the adaptor match
        my $length = length($processed_seq); # Get the length of the processed sequence
        $processed_count_3++; # We've matched the 3' adaptor so increment this counter
        if ($processed){ # If we have also matched the 5' adaptor (or no 5' adaptor was given) then:
          $processed_count_both++; # We matched and removed both adaptor sequences so the sequence we now have is fully processed and we can increment this count
          if (($length >= $minsize) && ($length <= $maxsize)){ # If it falls into the specified size class (default 18-30nt) then:
            print FASTA_FINAL ">$id\n$processed_seq\n"; # print the processed sequence out
            $processed_count_valid++; # Increment this counter
            $lengths{$length}++; # Increment a counter for the specified size class
          }
          next; # We've already found a match to both adaptors in the reverse complement sequence - therefore we can skip looking in the sense orientation
        }
      }
    }
    $processed = 0; # Reset this value 
    if (!$adaptor_5){ # If user not specified a 5' adaptor then increment
    $processed++;
    }
    elsif ($sequence =~m/$adaptor_5(\w+)$/){ # otherwise if we match the 5' adaptor 
    $sequence = $1; # Capture the sequence
    $processed_count_5++; # Increment
    $processed++; # Increment
    }
    if ($sequence =~m/^(\w+)($adaptor_3\w*)$/){ # Match 3' adaptor
    my $processed_seq = $1; # Capture sequence
    my $adaptor = $2;
    my $length = length($processed_seq); # Get length
    $processed_count_3++; # Increment counter
      if ($processed){ # If we've already removed 5' adaptor (or not given)
        $processed_count_both++; # We've successfully matched both adaptors
        if (($length >= $minsize) && ($length <= $maxsize)){ # Does it fall within specified size class?
          print FASTA_FINAL ">$id\n$processed_seq\n"; # if so print it out
          $processed_count_valid++; # increment counter
          $lengths{$length}++; # keep track of sRNA size counts
        }
      }
  }
}

close FASTA_FINAL;

}; # eval
if ($@) {
  SrnaTools::Exception::ModuleExecution->throw("$module_name, error while removing adaptors: $@") ;
}

open (OUTPUT, ">$output_dir/sequence_info.txt") or SrnaTools::Exception::ModuleExecution->throw("$module_name, could not open output file");

print OUTPUT "Sequences in input FASTQ/FASTA file: $total_count\n"; # Total number of input sequences
if (!$adaptor_5){ # If only 3' adaptor specified then only give these stats
  print OUTPUT "Sequences matching 3' adaptor sequence: $processed_count_3\n";
  print OUTPUT "Sequences matching 3' adaptor sequences within size range specified: $processed_count_valid\n\n";
}
else{ # Else give stats for 3' and 5' adaptor matches
  print OUTPUT "Sequences matching 5' adaptor sequence: $processed_count_5\n";
  print OUTPUT "Sequences matching 3' adaptor sequence: $processed_count_3\n";
  print OUTPUT "Sequences matching 3' and 5' adaptor sequences: $processed_count_both\n";
  print OUTPUT "Sequences matching 3' and 5' adaptor sequences within size range specified: $processed_count_valid\n\n";
}
print OUTPUT "sRNA length distribution (within specified size range)\n\n"; 

foreach my $element (sort keys %lengths){ # Outputs a size distribution of sRNAs that are valid (match adaptors specified and within size range)
  my $count = $lengths{$element};
  print OUTPUT "$element\t$count\n";
}
close OUTPUT;

return 1;
} #run


#########################################
# This replaces the fastq2fasta.pl script
#########################################
sub convert_fastq_to_fasta{
  my $self = shift ;
  my $infile = shift or die "need an input file name";
  my $outfile = $infile.'.converted';
  
  my $fasta=0;
  open (INPUT, "$infile") or die "Can't open input file\n" ;
  open (OUTPUT, ">", $outfile) or die "Can't open output file\n";
  my $id; # The sequence ID
  my $nextseq = 0;
  my $firstseq = 0;
  while (<INPUT>){
    if ($fasta){
      print OUTPUT $_;
    }
    elsif ($_=~m/^\@(\S+)/){ # If the line starts with an "@" symbol then we have a FASTQ file and any non-whitespace after the @ symbol must be the sequence ID which we capture - we then expect the sequence to be on the following line
    $id = $1; # capture the sequence ID
    chomp $id;
    $nextseq=1;
    $firstseq++;
    }
    elsif (($_=~m/^\>(\S+)/) && ($nextseq == 0) && ($firstseq == 0)){ # If however the line starts with a ">" symbol then the file must be in FASTA format and the following line will be the sequence
    print OUTPUT $_;
    $fasta++;
    $firstseq++;
    }
    elsif ($nextseq == 1){ # If the previous line began with either an "@" or a ">" then we know that this line should be the sequence
    my $sequence = $_; # capture the sequence
    chomp $sequence;
    $sequence =~ tr/a-z/A-Z/; # Make sequence uppercase
    $sequence =~ s/U/T/g; # Convert any 'U's to 'T's (RNA -> DNA)
    print OUTPUT ">$id\n$sequence\n"; # print out FASTA format: ">ID\nSequence\n"
    $nextseq=0; # This line was not the ID line so we know that the next line will not be the sequence line
    }
  }
  close INPUT;
  close OUTPUT;
  eval { unlink $infile } ;
  return $outfile ;
} # convert_fastq_to_fasta

1;