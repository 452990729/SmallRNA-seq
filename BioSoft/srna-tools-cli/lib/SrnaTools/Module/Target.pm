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


####################################################################
#
# This is a module of the sRNA tool kit.
# It finds targets for miRNAs using TXT followed by filtering 
# according to known target interaction criteria
#
####################################################################

package SrnaTools::Module::Target ;
use base SrnaTools::Module ;
use strict ;
use warnings;
use Cwd;

sub run{
  my $self = shift;

  require Bio::SeqIO ;
  require Bio::SearchIO;
  
  my $cwd = getcwd;
  
  my $module_name = "Target-finder";

  my $working_dir =  $self->job->job_working_dir;
  my $transcriptome_file_name = $self->param('transcriptome_file');
  my $transcriptome_name = $self->param('transcriptome_name');
  my $public_data_dir = $self->job->app->path_to('data_dir',$self->job->execution_env);
  my $transcriptome_dir = $public_data_dir.'/'. $self->job->transcriptomes_subfolder ;
  my $transcriptome_file = $transcriptome_dir.'/'.$transcriptome_file_name;
  my $data_dir = $working_dir.'/data/';
  my $results_dir = $working_dir.'/results/';
  my $infile = $data_dir.'/'.$self->param('infile');
  
  if (! -d $working_dir) {
    SrnaTools::Exception::ModuleExecution->throw(message=>"Could not access job working directory",log_msg=>"trying to access: $working_dir") ;
  }
  if (! -f $transcriptome_file || ! -r $transcriptome_file){
    SrnaTools::Exception::ModuleExecution->throw(message=>"Could not read transcriptome file.",log_msg=>"trying to read: $transcriptome_file") ;
  }
  if (! -r $infile){
    SrnaTools::Exception::ModuleExecution->throw(message=>"Could not read sRNA file.",log_msg=>"trying to read: $infile") ;
  }
  
  # check binaries
  SrnaTools::Exception::ModuleExecution->throw("fasta bin file not accessible or not executable") unless $self->binary_in_path('fasta34');
  
  SrnaTools::Exception::ModuleExecution->throw("clustalw2 bin file not accessible or not executable") unless $self->binary_in_path('clustalw2');
  
  SrnaTools::Exception::ModuleExecution->throw("RNAcofold bin file not accessible or not executable") unless $self->binary_in_path('RNAcofold');
  
  open (TXT, ">", "$results_dir/results.txt") or SrnaTools::Exception::ModuleExecution->throw("Can't create .fasta file\n");
  
  open (CSV, ">$results_dir/results.csv") or  SrnaTools::Exception::ModuleExecution->throw("Can't create .csv file\n");
  
  print CSV "sRNA ID, target gene accession, start-end position of target, target description\n";

  my %seqs;
  my %db_seqs;
  my %descs;
  my %fullacc;
  my @accs;
  
  $self->update_status("$module_name: running target identifcation") ;
  my $inseq = Bio::SeqIO->new('-file' => "<$infile",
                            '-format' => 'fasta' ) ;
  while (my $seq_obj = $inseq->next_seq ) {
  my $id = $seq_obj->id ;
  my $seq = $seq_obj->seq ;
  $seq=~ s/U/T/g;
  chomp ($id, $seq);
  $seqs{$id} = $seq;
  push (@accs, $id);
  }


  my $dbseq = Bio::SeqIO->new('-file' => "$transcriptome_file", '-format' => 'fasta' ) ;
  while (my $seq_obj = $dbseq->next_seq ) {
    my $id = $seq_obj->id ;
    my $seq = $seq_obj->seq ;
    my $desc = $seq_obj->description();
    chomp ($id, $seq, $desc);
    if ($id =~ m/(.*)\.\d+$/) {
      $db_seqs{$1} = $seq;
      $descs{$1}=$desc;
      $fullacc{$1} = $id;
    } else{
      $db_seqs{$id} = $seq;
      $descs{$id}=$desc;
      $fullacc{$id} = $id;
    }
  }

  my $number_of_seqs =  keys %seqs; 
  my $count = 0;

  foreach my $element (@accs) {    
    my $seq = $seqs{$element};
    my $id = $element;
    $count++;    
    
    
    
    my $tmp_seq_file = "$data_dir/tmpseq.$$";
    open (FILE, ">$tmp_seq_file") or SrnaTools::Exception::ModuleExecution->throw("Failed to open file: $tmp_seq_file");
    print FILE ">$seq\n$seq\n";
    close FILE;
    
    my $outfile = "$data_dir/tmp.$$";
    open (FILE, ">", $outfile) or SrnaTools::Exception::ModuleExecution->throw("Failed to create file: $outfile");
    close FILE;
    
    my $cmd = "fasta34 -W 3 -E 100 -b 500 -n -O '$outfile' '$tmp_seq_file' '$transcriptome_file'";
    my $execute = `$cmd`;
    my $targets_found = 0;
    
    my $parser = new Bio::SearchIO(-format => 'fasta', -file => "$outfile") ;
    while (my $result = $parser->next_result) {
      while(my $hit = $result->next_hit) {
        while (my $hsp = $hit->next_hsp) {
          my $hit_start = $hsp->start('hit') ;
          my $hit_end = $hsp->end('hit') ;
          my $query_start = $hsp->start('query') ;
          my $query_end = $hsp->end('query') ;
          my $strand = $hsp->strand('query') *  $hsp->strand('hit');
          my $query_string = $hsp->query_string;
          my $hit_string = $hsp->hit_string ;
          my $frac_identical = $hsp->frac_identical;
          my $accession = $hit->accession;
          my $description = $descs{$accession};
          my $est_seq = $db_seqs{$accession};
          my $est_length = length($est_seq);
          my $substr_start = $hit_start-15;
          
          if ($substr_start >= $est_length){
            next;
          }
          
          if ($substr_start < 0){
            $substr_start = 0;
          }
          
          my $end = 40;
          if (($substr_start+$end) >= $est_length){
            $end = $est_length-$substr_start;
          }
          
          my $est_subseq = uc(substr($est_seq, $substr_start, $end));
          my $est_subseq_rc = revcomp($est_subseq);
          if (length($est_subseq) < 14){
            next;
          }
          
          # Got to cd into data dir because clustalw2 can't
          # handle abs path 
          my $cwd = getcwd;
          chdir "$data_dir" or SrnaTools::Exception::ModuleExecution->throw("Failed to cd into $data_dir");
    
          # Make an alignment of RNA sequence vs EST target
          my $temp_align_file = "align_tmp";
          open (TMP, ">",$temp_align_file) or SrnaTools::Exception::ModuleExecution->throw("Failed to open $temp_align_file in ". `pwd`);
          
          print TMP ">RNA\n$seq\n";
          if ($strand == 1) {
            print TMP ">$accession/$hit_start-$hit_end\n$est_subseq\n";
          } elsif ($strand == -1) {
            print TMP ">$accession/$hit_start-$hit_end\n$est_subseq_rc\n";
          }
          close TMP;
          
          my $clustal =  `clustalw2 $temp_align_file`;
          open (TMP, "<","$temp_align_file.aln") or SrnaTools::Exception::ModuleExecution->throw("Failed to open $temp_align_file.aln file - clustalw may not have returned any resutls");
          
          my ($target_aln_string, $RNA_aln_string, $consensus_aln_string) = ("", "", "") ;
          while (<TMP>) {
            chomp;
            if (m/(RNA\s+)([A-Z\-]+)/){
              $RNA_aln_string = rev($2);
              $RNA_aln_string =~ m/^(\-*)([A-Z\-]+?)(\-*)$/;
              my ($lead, $seq, $tail) = ($1, $2, $3);
              $lead =~ s/\-/ /g;
              $tail =~ s/\-/ /g;
              $RNA_aln_string = $lead.$seq.$tail;
            } elsif ( m/(\S+\s+)([A-Z\-]+)/){
              $target_aln_string = revcomp($2);
            }
          }
          
          close TMP;
          unlink "$temp_align_file.aln";
          unlink "$temp_align_file";
          unlink "$temp_align_file.dnd";
          unlink "$tmp_seq_file";
          unlink "$outfile";
          
          chdir $cwd; # back to where we came from
          
          ### Calculate MFEs
          my ($heteroduplex_MFE,$homoduplex_MFE) = $self->calculate_MFE($target_aln_string,$RNA_aln_string); 
          ### Generate an 'alignment'
          foreach my $i (0 .. (length($target_aln_string)-1)){
            if (comp(substr($target_aln_string, $i, 1)) eq comp(comp(substr($RNA_aln_string, $i, 1)))) {
              $consensus_aln_string .= "|";
            }
            elsif (substr($target_aln_string, $i, 1) eq "G" and substr($RNA_aln_string, $i, 1) eq "T" ) {
              $consensus_aln_string .= "o";
            }
            elsif (substr($target_aln_string, $i, 1) eq "T" and substr($RNA_aln_string, $i, 1) eq "G" ) {
              $consensus_aln_string .= "o";
            }
            elsif (substr($RNA_aln_string, $i, 1) eq " " ) {
              $consensus_aln_string .= ".";
            }
            else {
              $consensus_aln_string .= " ";
            }
          }
          $consensus_aln_string = ".$consensus_aln_string.";

          ### Convert Ts to Us (DNA->RNA)
          $RNA_aln_string =~ s/T/U/g;
          $target_aln_string =~ s/T/U/g;
                            
          ### Test whether this predicted target satisfies the criteria
          my $target_ok = 1;
              
          ### Total number of mismatches must not exceed 4; G-U counts a 0.5 mismatch
          my @n_mismatches = ( $consensus_aln_string =~ m/\-|\s/g );
          my @n_GU = ( $consensus_aln_string =~ m/o/g );
          my $n_GU = @n_GU;
          my $mismatches = @n_mismatches + (0.5 * @n_GU);
          if ($mismatches > 4 ) {
            $target_ok = 0;
          }
      
          ### No more than one bulge in target (i.e. one gap in RNA)
          my @bulges = ($RNA_aln_string =~ /\-/g);
                  if (@bulges > 1 ) {
                  $target_ok = 0;
                  }
              
          ### No more than one bulge in RNA  (i.e. no gaps in target)
          @bulges = ($target_aln_string =~ /\-/g);
          if (@bulges > 1) {
            $target_ok = 0;
          }
              
          ### No more than two adjacent mismatches allowed
          if ($consensus_aln_string =~ m/[\so][\so][\so]/) {
            $target_ok = 0;
          }
                            
          ### No two adjacent mismatches in positions 2-12 (5') of RNA
          $consensus_aln_string =~ m/\.+(.+?)\.+/;
          my $subregion = substr($1, -12, 11);
          if ($subregion =~ /  |o | o/) {
            $target_ok = 0;
          }
              
          ### No mismatch in position 10 and 11
          $consensus_aln_string =~ m/\.+(.+?)\.+/;
          $subregion = substr($1, -11, 2);
          if ($subregion =~ / |o/) {
            $target_ok = 0;
          }
              
          ### No more than 2.5 mismatches in position 1-12
          $consensus_aln_string =~ m/\.+(.+?)\.+/;
          $subregion = substr($1, -12, 12);
          my @n_mismatches_in_this_subregion = ( $subregion =~ m/\-|\s/g );
          my @n_GU_in_this_subregion = ( $subregion =~ m/o/g );
          my $mismatches_in_this_subregion = @n_mismatches_in_this_subregion + 0.5 * @n_GU_in_this_subregion;
          if ($mismatches_in_this_subregion > 2.5) {
            $target_ok = 0;
          }

          $consensus_aln_string =~ s/\./ /g;
          my $alignment = "5' $target_aln_string 3'\n  $consensus_aln_string\n3' $RNA_aln_string 5'\n";
          my $satisfies_criteria = "";
          my $ratio = 0;
          
          unless (!$heteroduplex_MFE){
            $ratio = ($heteroduplex_MFE/$homoduplex_MFE);
          }
          
          if (( $target_ok ) && ($ratio >= 0.73)){
            my $estseq = $db_seqs{$accession};
            my $ori;
            if ($strand == -1){
              $ori = "Correct orientation";
            } else{
              $ori = "Wrong orientation";
              next;
            }
             
            my $accession = $fullacc{$accession};
            print TXT ">$id\t$accession\/$hit_start\-$hit_end\t$description\n$alignment\n";
            print  TXT "\n>$accession\n$estseq\n\n";
            print CSV "$id,$accession,$hit_start\-$hit_end,\"$description\"\n";
            $targets_found++;
          } else {
            next;
          }
        }
      }
    }
    if (!$targets_found){
      print TXT "No target found for sequence $id\n\n";
      print CSV "$id,No target found\n";
    }
  }

  close CSV;
  close TXT;
  chdir $cwd; # back to where we came from
  return 1;
}

sub comp {
  my $string = shift @_ or SrnaTools::Exception::ModuleExecution->throw("Called comp on null string!") ;
  my $seq = "" ;
  my @array = split(//, $string) ;
  while (my $char = shift @array) {
    if ($string =~ /U/ ) {
      ### RNA
      $char =~ tr/ACGU/UGCA/;
    } else {
      ### DNA
      $char =~ tr/ACGTU/TGCAA/;
    }
    $seq = $seq.$char 
  }
  return $seq; 
}

sub rev {
  my $string = shift @_ or SrnaTools::Exception::ModuleExecution->throw("Called rev on null string!") ;
  my $seq = "" ;
  my @array = split(//, $string) ;
  while (my $char = pop @array) {
    $seq = $seq.$char 
  }
  return $seq 
}

sub revcomp {
  my $string = shift or SrnaTools::Exception::ModuleExecution->throw("Called revcomp on null string!") ;
  my $seq = "" ;
  my @array = split(//, $string) ;
  while (my $char = pop @array) {
    if ($string =~ /U/ ) {
    ### RNA
    $char =~ tr/ACGU/UGCA/;
    } else {
      ### DNA
      $char =~ tr/ACGTU/TGCAA/;
    }
    $seq = $seq.$char;
  }
  return $seq;
}

sub calculate_MFE {
  my $self = shift;
  my $target_aln_string = shift or SrnaTools::Exception::ModuleExecution->throw("You did not specify a target aln string");
  my $RNA_aln_string = shift or SrnaTools::Exception::ModuleExecution->throw("You did not specify a RNA aln string");
  $target_aln_string =~ s/T/U/g;
  $RNA_aln_string =~ s/T/U/g;
  my ($heteroduplex_MFE, $homoduplex_MFE);
  ### Remove the dangling ends
  $RNA_aln_string =~ m/^(\s*)([ACGTUN\-]+?)(\s*)$/;
  my ($lead, $rna_seq_trimmed, $tail) = ($1, $2, $3);
  my $target_seq_trimmed = substr($target_aln_string, length($lead), length($rna_seq_trimmed));
  ### Connect the RNA and target with a linker sequence
  my $heteroduplex = rev($rna_seq_trimmed)."&".$target_seq_trimmed;
  my $homoduplex = rev($rna_seq_trimmed)."&".comp($rna_seq_trimmed);
  $heteroduplex =~ s/-//g;
  $homoduplex =~ s/-//g;
  
  my $working_dir =  $self->job->job_working_dir;
  my $data_dir = $working_dir.'/data/';
  my $tmp_file = "$data_dir/rna.tmp";
  
  open (TMP, ">$tmp_file") or SrnaTools::Exception::ModuleExecution->throw("Failed to open file: $tmp_file");
  
  print TMP ">heteroduplex\n$heteroduplex";
  close TMP;
  
  my $rna_cmd = "RNAcofold < $tmp_file";
  
  my $rna_execute = `$rna_cmd`;
  if ($rna_execute =~ /\((\s*\-[\d\.]+)\)/) {
    $heteroduplex_MFE = $1;
  }
  open (TMP, ">$tmp_file") or SrnaTools::Exception::ModuleExecution->throw("Failed to open file: $tmp_file");
  
  print TMP ">homoduplex\n$homoduplex";
  close TMP;
  $rna_cmd = "RNAcofold < $tmp_file";
  $rna_execute = `$rna_cmd`;
  if ($rna_execute =~ /\((\s*\-[\d\.]+)\)/) {
    $homoduplex_MFE = $1;
  }
  unlink "$tmp_file";
  unlink "heteroduplex_ss.ps";
  unlink "homooduplex_ss.ps";
  return ($heteroduplex_MFE,$homoduplex_MFE);
}


1;