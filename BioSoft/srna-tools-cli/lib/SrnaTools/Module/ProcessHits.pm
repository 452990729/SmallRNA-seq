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
# This is a module of the sRNA tool kit and part of the Mircat 
# pipeline. This class does the actual miRNA idenitfication.
#
# TODO This class has been converted from a standalone script and might
# still need some cleanup (e.g. passing of parameters)
#
#####################################################################

package SrnaTools::Module::ProcessHits ;
use base SrnaTools::Module ;
use strict;
use warnings;

my (%chr_length);

sub run{
my $self = shift ;

# modules that might be in non-standard
# locations must be dealt with here
# See SrnaTools::Module for details
require Bio::Seq::SeqFactory;
require Bio::Seq;
require Bio::DB::Fasta;
my $randfold =''; # Please enter the full path to randfold in the quotes unless it's in your path - obtain from http://bioinformatics.psb.ugent.be/software.php
my $rnafold = ''; # Please enter the full path to RNAfold in the quotes unless it's in your path - obtain from http://www.tbi.univie.ac.at/~ivo/RNA/

my $module_name = "miRCat" ;
my $mirs_found = 0;
my (%seen, %seen_acc, %abundance);
my ($lib_dir, $dir, $input, $no_complex_loops, $min_hairpin_len, $clust_length, $create_gff, $max_overlap_length, $genome, $help, $ext, $min_mirna_hits, $withrandfold, $min_paired, $min_hits, $max_gaps, $max_genome_hits, $min_mature, $max_mature, $min_gc, $max_unpaired, $ori_filter);

$withrandfold = 1;

my $minenergy = $self->param('min_energy') ;
$genome = $self->param('genome_file');
$ext = $self->param('window_length');
$min_paired  = $self->param('min_paired');
$min_hits = $self->param('min_abundance');
$max_gaps = $self->param('max_gaps');
$max_genome_hits = $self->param('genomehits');
$min_mature = $self->param('minsize');
$max_mature = $self->param('maxsize');
$min_gc = $self->param('min_gc') ;
$max_unpaired = $self->param('max_percent_unpaired') ;
$max_overlap_length = $self->param('max_overlap_length');
$clust_length = $self->param('max_unique_hits');
$ori_filter = $self->param('percent_orientation');
$min_hairpin_len = $self->param('min_hairpin_len') ;
$no_complex_loops =  $self->param('no_complex_loops'); 
my $hit_dist =  $self->param('hit_dist')||200;
my $pval =  $self->param('pval');
$dir =  $self->job->job_working_dir ;
$lib_dir = $self->job->app->path_to('lib_dir',$self->job->execution_env); 

if( not $genome or not $dir ) {
  SrnaTools::Exception::ModuleExecution->throw(message=>"$module_name, no genome or directory given");
}

if (!$randfold){
$randfold = "randfold";
SrnaTools::Exception::ModuleExecution->throw("randfold bin file not accessible or not executable") unless $self->binary_in_path($randfold);
}

if (!$rnafold){
$rnafold = "RNAfold";
SrnaTools::Exception::ModuleExecution->throw("RNAfold bin file not accessible or not executable") unless $self->binary_in_path($rnafold);
}
if ($minenergy > 0){
$minenergy = -25.0;
}
if (!$min_hairpin_len){
$min_hairpin_len = 75;
}


if (!$ori_filter){
$ori_filter = 0;
}
chomp ($genome, $dir, $lib_dir);
if (system("$lib_dir/create_fasta_db.pl $genome")!=0 ) {
  SrnaTools::Exception::ModuleExecution->throw(message=>"$module_name, Fasta DB creation step failed",log_msg=>"tried to run:$lib_dir/create_fasta_db.pl $genome ");
}
my $seqinx = Bio::DB::Fasta->new("$genome");

###############################
# Die if something is wrong.. #
###############################

if (!-e $genome){
  SrnaTools::Exception::ModuleExecution->throw(message=>"$module_name, couldn't locate genome file",log_msg=>"file: $genome");
}


################################################
# Use default values if they aren't specified: #
################################################
# TODO this should not be needed anymore because the App passes defaults
if (!$ext){
warn "Extention length not specified - using default value of 100nt\n";
$ext = 150;
}
if (!$min_paired){
warn "Minimum number of bases in miRNA to be paired not specified - using default value of 18\n";
$min_paired = 17;
}
if (!$min_hits){
warn "Minimum number of hits needed to test for miRNA not given - using default value of 1\n";
$min_hits = 1;
}
if (!$max_gaps){
warn "Maximum number of consecutive mismatches allowed between candidate miRNA and miRNA* - using default value of 3\n";
$max_gaps = 3;
}
if (!$max_genome_hits){
warn "Maximum number of genomic loci of candidate mature miRNA not given - using default value of 16\n";
$max_genome_hits = 16;
}
if ((!$max_mature) || (!$min_mature)){
warn "Minimum/maximum length of mature miRNA not specified using defaults (min 18nt max 25nt)\n";
$min_mature = 20;
$max_mature = 22;
}
if (!$min_gc){
warn "Minimum G/C percentage of sRNA not specified - using default (20%)\n";
$min_gc = 10;
}
if (!$max_unpaired){
warn "Maximum percentage unpaired nucleotides in hairpin not specified - using default (50%)\n";
$max_unpaired = 50;
}
if (!$max_overlap_length){
$max_overlap_length = 70;
}

if (!$clust_length){
$clust_length = 3;
}

# Print to file
open (CSV, '>', $dir.'/results/results/output.csv') or SrnaTools::Exception::ModuleExecution->throw("Could not generate output file") ;

chomp ($min_hits);
if ($withrandfold){
  print CSV "Chromosome,Start,End,Orientation,Abundance,Sequence,sRNA length,# Genomic Hits,Hairpin Length,Hairpin % G/C content,Minimum Free Energy,Adjusted MFE,Randfold p-value,miRNA*\n";
}
else{
  print CSV "Chromosome,Start,End,Orientation,Abundance,Sequence,sRNA length,# Genomic Hits,Hairpin Length,Hairpin % G/C content,Minimum Free Energy,Adjusted MFE,miRNA*\n";
}

my %all_coords;
my %chromosomes;
my %genome_hits;
my $mature_length;
my %unique;
#######################
# Create output files #
#######################
my $file = "output";

open (STRUCTURE, ">$dir/results/results/miRNA_hairpins.txt")  or SrnaTools::Exception::ModuleExecution->throw("Could not generate miRNA_hairpins.txt file") ;
if ($create_gff){
open (GFF, ">$dir/results/results/$file.gff")  or SrnaTools::Exception::ModuleExecution->throw("Could not generate gff file") ;
}
open (FASTA, ">$dir/results/results/FASTA.txt")  or SrnaTools::Exception::ModuleExecution->throw("Could not generate FASTA.txt file") ;
my $inputfile = "$dir/data/patman.out";
open (INPUT, "$inputfile")  or SrnaTools::Exception::ModuleExecution->throw(message=>"Could not open genomic coordinates file",log_msg=>"tried to open $inputfile") ;

################################################
# Read in genomic coordinates of all sRNA hits #
################################################

# Get non-redundant sequences
my %locs;
my %unqiue;
while (<INPUT>){
  if ($_=~m/^(\S+)\s*.*\t([AGCTNU]+?)\_\S+\((\d+)\)\t(\d+)\t(\d+)\t(\S+)\t\d+/i){
    my $chromosome = $1;
    my $start = $4;
    my $end = $5;
    my $strand = $6;
    my $sequence = $2;
    my $abundance = $3;
    chomp ($chromosome,$start,$end,$strand,$sequence,$abundance);
    $sequence =~ tr/U/T/;
    $sequence =~ tr/a-z/A-Z/;
    if (!$locs{"$chromosome\/$start\-$end$strand"}){
      $genome_hits{$sequence}++;
      $locs{"$chromosome\/$start\-$end$strand"}++;
    }
    if (!$chromosomes{$chromosome}){
      $chromosomes{$chromosome}++;
    }
      $abundance{$sequence} = $abundance;
      $unique{$sequence}++;
  }
}
undef %locs;
close INPUT;
my $count =0;

# Goes through hits chromosome by chromosome:
my %chr_coords;
foreach my $chr (sort keys %chromosomes){
  open (INPUT, "$inputfile") or SrnaTools::Exception::ModuleExecution->throw(message=>"Could not open genomic coordinates file",log_msg=>"tried to open $inputfile") ;

  my $scount=0;
  while (<INPUT>){
    if ($_=~m/^(\S+)\s*.*\t([AGCTNU]+?)\_\S+\((\d+)\)\t(\d+)\t(\d+)\t(\S+)\t\d+/){
      my $chromosome = $1;
      my $start = $4;
      my $end = $5;
      my $strand = $6;
      my $sequence = $2;
      my $abundance = $3;
      $sequence =~ tr/U/T/;
      $sequence =~ tr/a-z/A-Z/;
      chomp ($chromosome, $start, $end, $strand, $sequence);
      my $se = "$start\-$end$strand";
      if ($strand eq '+'){
        $strand = "1";
      } else{
        $strand = "-1";
      }
      my $no_hits = $abundance{$sequence};
      if ($chromosome eq $chr){
        $scount++;
        if (!$chr_coords{$se}){
          $chr_coords{$se}= {'start' =>$start, 'end' =>$end, 'strand' => $strand, 'chromosome' => $chromosome, 'no_hits' =>$no_hits, 'sequence' =>$sequence};
        }
      }
    } 
  } #while
  close INPUT;
  my $last_end = 0;
  my $cluster_length = 0;
  my $cluster_count = 0;
  my %used_start_positions;

  ######################
  # Check for clusters #
  ######################

  my $mixed_no = 1;
  my %mix_count;
  my $plus = 0;
  my $minus = 0;
  my %hitcount;
  my $lastseq;
  my $overlap_length = 0;
  my $max_overlap = 0;
  my %seen_seq_ori;
  
  foreach my $element (sort {$chr_coords{$a}->{'start'} <=> $chr_coords{$b}->{'start'}} keys %chr_coords){
    my $start = $chr_coords{$element}->{'start'};
    my $end = $chr_coords{$element}->{'end'};
    my $strand = $chr_coords{$element}->{'strand'};
    my $sequence = $chr_coords{$element}->{'sequence'};
    my $seq_length = length($sequence);
    my $no_hits = $chr_coords{$element}->{'no_hits'};
    my $o_id = $chr_coords{$element}->{'sequence'};
    my $dist = $start - $last_end;
    
    if ($last_end == 0){
      $hitcount{"$chr\t$start\-$end($no_hits)\t$strand\t$sequence $o_id"}=$no_hits;
      $cluster_count++;
      $overlap_length = $seq_length;
    } elsif (($last_end >= $start)){
      $hitcount{"$chr\t$start\-$end($no_hits)\t$strand\t$sequence $o_id"}=$no_hits;
      if ($end > $last_end){
        my $overlap = ($seq_length-($last_end-$start))-1;
        $overlap_length = $overlap_length+$overlap;
        if ($overlap_length > $max_overlap){
          $max_overlap = $overlap_length;
        }
      }
    } elsif ($dist <= $hit_dist){
      $hitcount{"$chr\t$start\-$end($no_hits)\t$strand\t$sequence $o_id"}=$no_hits;
      $cluster_count++;
      
      if ($overlap_length > $max_overlap){
        $max_overlap = $overlap_length;
      }
      $overlap_length = $seq_length;
    } elsif ($cluster_count <= $clust_length){
      my $maxhit = 0;
      my $maxstrand;
      my $maxseq;
      my $out ="";
      
      if ($overlap_length > $max_overlap){
        $max_overlap = $overlap_length;
      }
      $overlap_length = $seq_length;
      
      foreach my $details (sort keys %hitcount){
        chomp $details;
        my $element = $hitcount{$details};
        chomp $element;
        if ($details=~m/^\S+\s+(\d+)\-(\d+)\(\d+\)\s+(\-?\d+)\s+(\S+).+/){
          if ($element > $maxhit){
            my $startpos = $1;
            my $endpos = $2;
            my $strand = $3;
            my $seq_string = $4;
            chomp ($startpos,$endpos, $strand,$seq_string);
            $out = "$startpos$endpos";
            $maxhit = $element;
            $maxstrand = $strand;
            $maxseq = $seq_string;
          }
        }
      }
                
      $cluster_count = 0;
      my $strand_bias = 0;
      if ($maxstrand == 1){
        if ($minus == 0){
        $strand_bias = 100;
        }
        else{
        $strand_bias = ($plus/($minus+$plus))*100;
        }
      } else{
        if ($plus == 0){
        $strand_bias = 100;
        }
        else{
        $strand_bias = ($minus/($minus+$plus))*100;
        }	
      }
      
      if (($genome_hits{$maxseq} <= $max_genome_hits) && ($strand_bias >= $ori_filter) && ($max_overlap <= $max_overlap_length)){
        $used_start_positions{$out}=1;
      }
      $plus = 0;
      $minus = 0;
      undef %hitcount;
      undef %seen_seq_ori;
      $max_overlap = $seq_length;
      $hitcount{"$chr $start\-$end($chr_coords{$element}->{'no_hits'})\t$strand\t$sequence $o_id"}=$no_hits;
    } else{
                $cluster_count = 0;
		$plus = 0;
		$minus = 0;
		undef %hitcount;
		undef %seen_seq_ori;
		$overlap_length = $seq_length;
		$max_overlap = $seq_length;
		$hitcount{"$chr $start\-$end($chr_coords{$element}->{'no_hits'})\t$strand\t$sequence $o_id"}=$no_hits;
                $cluster_count++;
                }
		if ($strand == 1){
			if (!$seen_seq_ori{$sequence}){
			$plus = $plus+$no_hits;
			$seen_seq_ori{$sequence}++
			}
		}
		else{
			if (!$seen_seq_ori{$sequence}){
			$minus = $minus+$no_hits;
			$seen_seq_ori{$sequence}++
			}
		}
		unless ($end < $last_end){
		$last_end = $end;
		}
	$lastseq = $sequence;
        }
        if ($cluster_count <= $clust_length){
	my $maxhit = 0;
	my $maxstrand;
	my $maxseq;
	my $out ="";
		if ($overlap_length > $max_overlap){
		$max_overlap = $overlap_length;
		}
		foreach my $details (sort keys %hitcount){
		chomp $details;
		my $element = $hitcount{$details};
		chomp $element;
			if ($details=~m/^\S+\s+(\d+)\-(\d+)\(\d+\)\s+(\-?\d+)\s+(\S+).+/){
				if ($element > $maxhit){
                                my $startpos = $1;
				my $endpos = $2;
				my $strand = $3;
				my $seq_string = $4;
                                chomp ($startpos,$endpos, $strand,$seq_string);
				$out = "$startpos$endpos";
				$maxhit = $element;
				$maxstrand = $strand;
				$maxseq = $seq_string;
                                }
			}
                }
		my $strand_bias = 0;
			if ($maxstrand == 1){
				if ($minus == 0){
				$strand_bias = 100;
				}
				else{
				$strand_bias = ($plus/($minus+$plus))*100;
				}
			}
			else{
				if ($plus == 0){
				$strand_bias = 100;
				}
				else{
				$strand_bias = ($minus/($minus+$plus))*100;
				}	
			}
		if (($genome_hits{$maxseq} <= $max_genome_hits) && ($strand_bias >= $ori_filter) && ($max_overlap <= $max_overlap_length)){
		$used_start_positions{$out}=1;
		}
        }
        $cluster_count = 0;
	undef %hitcount;
	undef %seen_seq_ori;
	$plus = 0;
	$minus = 0;
        $last_end = 0;
	my $tmp=0;

	
	
	####################
	# Check for miRNAs #
	####################
	foreach my $element (sort {$chr_coords{$a}->{'start'} <=> $chr_coords{$b}->{'start'}} keys %chr_coords){
	my $start = $chr_coords{$element}->{'start'};
	my $end = $chr_coords{$element}->{'end'};
	my $strand = $chr_coords{$element}->{'strand'};
	my $sequence = $chr_coords{$element}->{'sequence'};
	my $seq_length = length($sequence);
	my $hits = $chr_coords{$element}->{'no_hits'};
	my $o_id = $chr_coords{$element}->{'sequence'};
	chomp ($start, $end, $strand, $sequence, $hits, $o_id);
	my $gc_percentage = &get_gc($sequence); # get multiple hits that are very AT rich - filter these out (based on input settings)
	if (($hits >= $min_hits) && ($seq_length >= $min_mature) && ($seq_length <= $max_mature) && ($used_start_positions{"$start$end"}) && ($gc_percentage >= $min_gc)){
          $mature_length = length($sequence);
          my $retry = 1;
          my $newext1;
          my $newext2;
          my %outputdata;
          my %outputstructure;
          my %extend_lengths;
          print STDERR "Testing: $chr\/$start\-$end\t$sequence\n";
			while ($retry < 15){

			###############################
			# Use several folding windows #
			###############################
			
				if ($retry == 1){
				$newext1 = $ext;
				$newext2 = $ext;
				}
				elsif ($retry == 2){
				$newext1 = (int($ext*0.2));
				$newext2 = (int($ext*1.8));
				}
				elsif ($retry == 3){
				$newext1 = (int($ext*1.8));
				$newext2 = (int($ext*0.2));
				}
				elsif ($retry == 4){
				$newext1 = (int($ext*.75));
				$newext2 = (int($ext*.75));
				}
				elsif ($retry == 5){
				$newext1 = (int($ext*0.2));
				$newext2 = (int($ext*0.8));
				}
				elsif ($retry == 6){
				$newext1 = (int($ext*0.8));
				$newext2 = (int($ext*0.2));
				}
				elsif ($retry == 7){
				$newext1 = (int($ext*0.9));
				$newext2 = (int($ext*0.1));
				}
				elsif ($retry == 8){
				$newext1 = (int($ext*0.1));
				$newext2 = (int($ext*0.9));
				}
				elsif ($retry == 9){
				$newext1 = (int($ext*0.5));
				$newext2 = (int($ext*0.5));
				}
				elsif ($retry == 10){
				$newext1 = (int($ext*0.2));
				$newext2 = (int($ext*0.7));
				}
				elsif ($retry == 11){
				$newext1 = (int($ext*0.7));
				$newext2 = (int($ext*0.2));
				}
				elsif ($retry == 12){
				$newext1 = (int($ext*0.6));
				$newext2 = (int($ext*0.4));
				}
				elsif ($retry == 13){
				$newext1 = (int($ext*0.6));
				$newext2 = (int($ext*0.4));
				}
				elsif ($retry == 14){
				my $i = 0;
					foreach my $candidate (sort { $a <=> $b } keys(%outputdata)){
						if ($i == 0){
						my $outval = $outputdata{$candidate};
						my $structure = $outputstructure{$candidate};
						print CSV "$outval";
							if (!$seen{$o_id}){
							print FASTA ">$o_id\n$sequence\n";
							$seen{$o_id}++;
							}
						print STRUCTURE "$structure";
							if ($create_gff){
							print GFF "$chr\t\.\tmiRNA\t$start\t$end\t$hits\t$strand\t\.\t## Predicted miRNA\n";
							}
						$i++;
						}
					}
				}
			$retry++;
			my @info = &get_seq_ext($chr, $start, $end, $strand, $newext1, $newext2, $seqinx); # Extends sequence from genomic coordinates
			my $seq = $info[0];
			$newext1 = $info[1];
			$newext2 = $info[2];
			chomp ($seq, $newext1, $newext2);
				if (exists $extend_lengths{"$newext1$newext2"}){
				next;
				}
			$extend_lengths{"$newext1$newext2"}++;
			my $out = $seq->seq;
			my @structure = &getstructure($out, $newext1, $newext2,$dir,$rnafold,$mature_length); # Folds the extended sequence
			my $raw_structure = $structure[0];
				if ($raw_structure eq 0){
				next;
				}
			my $raw_seq = $structure[1];
			my $mirna_pos = $structure[2];
			chomp ($raw_structure, $raw_seq, $mirna_pos);
			my @ret = &trim_structure($raw_structure, $raw_seq, $mirna_pos,$max_gaps, $min_paired); # Trims the structure back to a hairpin
			my $trimmed_structure = $ret[0];
			my $trimmed_seq = $ret[1];
			my $trimmed_mirna_pos = $ret[2];
				if ($trimmed_structure eq 0){
				next;
				}
			chomp ($trimmed_structure, $trimmed_seq, $trimmed_mirna_pos);
			my @analysis = &analyse_structure($trimmed_structure,$min_hairpin_len,$no_complex_loops,$max_unpaired); # Analyse to see whether it forms a valid hairpin
			my $is_mirna = $analysis[0];
			my $hairpin_length = $analysis[1];
			chomp ($is_mirna, $hairpin_length);
				if ($is_mirna == 2){
				next;
				}
			my @returned = &refold($is_mirna, $hairpin_length, $trimmed_seq, $trimmed_structure,$dir,$withrandfold,$randfold,$rnafold); # Refold using constraints in order to get MFE of trimmed hairpin
			my $mfe = $returned[0];
			my $isgood = $returned[1];
			my $return_seq = $returned[2];
			my $return_struc = $returned[3];
			my $return_prob = $returned[4];
			my $return_length = $returned[5];
			my $gc_content = &get_gc($return_seq);
			my $amfe = ($mfe/$return_length)*100;
				if ((($isgood == 1) || ($isgood == 3)) && ($amfe <= $minenergy)){ # If it's a valid hairpin then print it out
					if ($strand eq '1'){$strand = '+';}
					elsif ($strand eq '-1'){$strand = '-';}
				my @mir_star_ret = &find_mir_star($start, $end, $return_seq, $trimmed_mirna_pos, $sequence);
				my $mirstar_valid = &validate_mir_star($start, $end, $return_seq, $trimmed_mirna_pos,$max_gaps,$min_paired);
				my $mirstar_present = "";
				my $count = 0;
					foreach my $mirstar (@mir_star_ret){
					chomp $mirstar;
						if (exists $genome_hits{$mirstar}){
						my $mirstar_abundance = $abundance{$mirstar};
						$mirstar_present = "$mirstar_present$mirstar($mirstar_abundance) ";
							if ($count == 0){
							my $structure_markup = &markup_mirstar($mirstar,$trimmed_mirna_pos,$return_seq);
							$trimmed_mirna_pos = $structure_markup;
							}
						$count++;
						}
					}
					if (!$mirstar_present){
					$mirstar_present = "NO";
					}
					my $srna_length = length($sequence);
					my $no_genome_hits = $genome_hits{$sequence};
					if ($withrandfold){
 						if (($return_prob < $pval) && ($mirstar_valid)){
						$outputdata{$mfe}="$chr,$start,$end,$strand,$hits,$sequence,$srna_length,$no_genome_hits,$return_length,$gc_content,$mfe,$amfe,$return_prob,$mirstar_present\n";
						$outputstructure{$mfe}=">$o_id"."_$chr\/$start\-$end\n$return_seq\n$trimmed_mirna_pos\n\n";
						$mirs_found++;
						}
					}
					elsif ($mirstar_valid){
					$outputdata{$mfe}="$chr,$start,$end,$strand,$hits,$sequence,$srna_length,$no_genome_hits,$return_length,$gc_content,$mfe,$amfe,$mirstar_present\n";
					$outputstructure{$mfe}=">$o_id"."_$chr\/$start\-$end\n$return_seq\n$trimmed_mirna_pos\n\n";
					$mirs_found++;
					}
				}
			}
		}
	}
	if ($create_gff){
		foreach my $element (keys %chr_coords){
		my $start = $chr_coords{$element}->{'start'};
		my $end = $chr_coords{$element}->{'end'};
		my $strand = $chr_coords{$element}->{'strand'};
		my $chromosome = $chr_coords{$element}->{'chromosome'};
		my $acc = $chr_coords{$element}->{'sequence'};
		my $abundance = $chr_coords{$element}->{'no_hits'};
		chomp ($start, $end, $strand, $chromosome, $acc, $abundance);
		if ($strand eq '1'){$strand = '+';}
		elsif ($strand eq '-1'){$strand = '-';}
			if (!$seen_acc{$acc}){
			print GFF "$chromosome\t\.\tmisc_RNA\t$start\t$end\t$abundance\t$strand\t\.\t$acc\n";
			}
		}
	}
undef %used_start_positions;
undef %chr_coords;
}
if ($create_gff){
close GFF;
}
if (!$mirs_found){
print CSV "No candidate miRNAs found - please make sure that you have uploaded a redundant input sRNA dataset\n";
}
close CSV ;
return 1;
}

sub markup_mirstar{
my $mirstar = shift;
my $structure = shift;
my @structure = split(//,$structure);
my $sequence = shift;
my $return_structure = "";
my $char = "";
my $star_structure = "";
	if ($structure =~m/>/){
	$char = '{'; 
	}
	else{
	$char = '}'; 	
	}
	if ($sequence=~m/(.*)($mirstar)(.*)/){
	my $head = $1;
	my $match = $2;
	my $tail = $3;
	chomp ($head, $match, $tail);
	my $first = length($head);
	my $middle = length($match);
	$middle=$middle+$first;
	my $last = length($tail);
	$last = $last+$middle;
	my $i = 0;
		foreach my $element (@structure){
		chomp $element;
			if ($i < $first){
			$return_structure = "$return_structure"."$element";
			}
			elsif ($i < $middle){
				if ($element eq '.'){
				$return_structure = "$return_structure"."=";
				$star_structure = "$star_structure"."=";
				}
				elsif (($element eq '(') || ($element eq ')')){
				$return_structure = "$return_structure"."$char";
				$star_structure = "$star_structure"."=";
				}
				else{
				$return_structure = "$return_structure"."$element";
				$star_structure = "$star_structure$element";
				}
			}
			else{
			$return_structure = "$return_structure"."$element";
			}
		$i++;
		}
	}
	if ($star_structure =~m/[\-\<\>]/){
	return($structure);
	}
	else{
	return($return_structure);
	}
}

# Extend sequence
sub get_seq_ext {
    my $id     = shift;
    my $start  = shift;
    my $end    = shift;
    my $strand = shift;
    my $extend1 = shift;
    my $extend2 = shift;
    my $seqinx = shift ;
    my ($header, $sequence, $s);
    my $seq = new Bio::Seq;
            unless (exists $chr_length{$id}){
                    eval {
                       $seq = $seqinx -> get_Seq_by_acc( $id );
                       my $clength = $seq->length();
                       $chr_length{$id} = $clength;
                    };
                    if( not $seq or $@ ) {
                    SrnaTools::Exception::ModuleExecution->throw("Sequence not found in seq DB" );
                           }
                   }
    my $oldstart = $start;
    $start = $start - $extend1;
    $end = $end + $extend2;
    
    my $length = $chr_length{$id};
            if( $start < 1 ) {
                $extend1 = $oldstart; # I think that this is wrong.. oldstart -1?
        $start = 1;
            }
            if( $end > $length ) {
                $extend2 = ($extend2 - ($end - $length));
        $end = $length;
            }

    my $newseq;
    my $truncseq = $seqinx -> subseq( $id, $start, $end );
    $newseq = Bio::Seq->new( -display_id => "$id\/$start-$end",
                                        -seq => $truncseq);
        if (($strand eq '-1')||($strand eq '-')){
                $newseq = $newseq -> revcom();
                my $tmp1 = $extend1;
                my $tmp2 = $extend2;
                $extend1 = $tmp2;
                $extend2 = $tmp1;
        }
my @return;
push (@return, $newseq, $extend1, $extend2);
return (@return);
}


# Takes extended sequence and folds
sub getstructure{
my $out = shift;
my $extend1=shift;
my $extend2=shift;
my $dir =shift;
my $rnafold = shift;
my $mature_length = shift;
chomp ($out, $extend1, $extend2);
my $isgood;
my $return_seq;
my $return_struc;
my $return_mirna_pos;
my $mirna_pos;
my $structure;
my $mfe = 0;
my @return;
open (FILE1, ">$dir/data/file1.$$") or SrnaTools::Exception::ModuleExecution->throw("Couldn't create temporary file in folding stage") ;
print FILE1 "$out";
close FILE1;

open (STRUC1, "$rnafold < $dir/data/file1.$$ |") or SrnaTools::Exception::ModuleExecution->throw("Couldn't perform RNAfold step" ) ;
	while (<STRUC1>){
		if ($_=~m/^([\.\(\)\<\>]+)\s+\((\s?-?\d*.\d+)\)$/){
		$mfe = $2;
		chomp $mfe;
		$structure = $1;
		 	unless ($structure =~m/[\(\)\<\>]+/){
			$structure =0;
			push (@return, $structure, $out, $mirna_pos);
			return (@return);
			}
		my @tmparry = split (//, $structure);
		my $p = 0;
			while ($p < $mature_length){
				if ($tmparry[$extend1+$p] eq '('){ 
				$tmparry[$extend1+$p]='<';
				}
				elsif ($tmparry[$extend1+$p] eq ')'){
				$tmparry[$extend1+$p]='>';
				}
				elsif ($tmparry[$extend1+$p] eq '.'){
				$tmparry[$extend1+$p]='-';
				}	
			$p++;
			}
			foreach my $el (@tmparry){
			$mirna_pos = "$mirna_pos$el";
			}
		}
		elsif ($_=~m/^\w+/){
		next;
		}
		else{
		$structure =0;
		$mirna_pos=0;
		}
		
	}
chomp ($structure, $mirna_pos);
push (@return, $structure);
push (@return, $out);
push (@return, $mirna_pos);
close STRUC1;
unlink "out.$$";
unlink "file1.$$";
return (@return);
}


# Refolds sequence
sub refold{
my $isgood = shift;
my $struc_length = shift;
my $return_seq = shift;
my $return_struc = shift;
my $dir = shift ;
my $withrandfold = shift;
my $randfold = shift;
my $rnafold = shift ;
my ($mfe, $prob) = 0;
my @return;
	if (($isgood == 1) || ($isgood == 3)){ # Need to get actual MFE again as current MFE is for whole sequence
	open (FILE2, ">$dir/data/file2.$$") or SrnaTools::Exception::ModuleExecution->throw("Couldn't create temp file in refolding stage" ) ;
	print FILE2 ">tmp\n";
	print FILE2 "$return_seq\n";
	close FILE2;
		if ($withrandfold){
			eval{
			system ("$randfold -d $dir/data/file2.$$ 100 > $dir/data/randfold.$$") == 0 or die; # Run randfold on all supposed good miRNAs and get a p-value
			};
			if($@) {
                          SrnaTools::Exception::ModuleExecution->throw("randfold step failed") ;
			}  
		unlink "$dir/data/file2.$$";
		open (RANDFOLD, "$dir/data/randfold.$$") or SrnaTools::Exception::ModuleExecution->throw("Couldn't open randfold results" ) ;
			while (<RANDFOLD>){
				if ($_=~m/^(\S+)\s+(\S+)\s+(\S+).*/){
				$mfe = $2;
				$prob = $3;
				chomp ($mfe, $prob);
				}
			}
		unlink "$dir/data/randfold.$$";
		}
		else{
		$prob = "N/A";
		open (FILE2, ">$dir/data/file2.$$") or SrnaTools::Exception::ModuleExecution->throw("Couldn't create temp file in refolding stage 2") ;
		print FILE2 "$return_seq\n";
		print FILE2 "$return_struc\n";
		close FILE2;
		open (STRUC2, "$rnafold -C < $dir/data/file2.$$ |") or SrnaTools::Exception::ModuleExecution->throw("Couldn't open refolded RNAfold output") ;;
			while (<STRUC2>){
				if ($_=~m/^([\.\(\)\<\>-]+)\s+\((\s?-?\d*.\d+)\)$/){
				$mfe = $2;
				chomp $mfe;
				}
			}
		close STRUC2;
		unlink "$dir/data/file2.$$";
		}
	}
unlink "$dir/data/file1.$$";
push (@return, $mfe);
push (@return, $isgood);
push (@return, $return_seq);
push (@return, $return_struc);
push (@return, $prob);
push (@return, $struc_length);
return (@return);
}


# Takes folded and extended sequence and trims down structure to try to form a hairpin
sub trim_structure{
my $structure = shift;
chomp $structure;
my $orig_seq = shift;
my $mirna_pos = shift;
my $max_gaps = shift ;
my $min_paired = shift ;
my @seq = split (//, $orig_seq);
my @struc = split (//, $structure);
my @mirnaposition = split (//, $mirna_pos);
my @pairs = @struc;
my $tmpsize = @pairs;
my $tmpsize2 = @mirnaposition;
my %opens;
my %closes;
my $last = '(';
my $i=0;
my $j = @struc;
my $countopen =0;
my @closing;
my $change =0;
my ($seq_start, $seq_end) = 0;
		# Find which columns pair with which
		while ($i < $j){
			if ($pairs[$i] eq '('){
			$countopen++;
			unshift (@closing, $countopen);
			$pairs[$i]=$countopen;
			$last = '(';
			}
			if ($pairs[$i] eq ')'){
			my $countclose = shift @closing;
			$pairs[$i]=$countclose;
			$countclose--;
			$last = ')';
			}
		$i++;
		}
	my $bases;
	my $counter = 0;
	# Mark up potential miRNA position on sequence 
	foreach my $el (@mirnaposition){
		if (($el eq '<')||($el eq '>')||($el eq '-')){
		$bases = "$bases$struc[$counter]";
		}
	$counter++;
	}
	$counter =0;
	chomp $bases;
	my $mir_length = length($bases);
	my $max_flank = 25-$mir_length;
	my $left_flank = $max_flank;
	my $right_flank = 0;
	my @mir_struc;


	# Get all 25nt windows containing potential mature miRNA
	$counter = 0;
	my @structure_lines;
	my @sequence_lines;
		if ($max_flank <= 0){
		push (@structure_lines, $bases);
		}
	my $found =0;
	my $startpos=0;
	while ($max_flank > 0){
		foreach my $el (@mirnaposition){
			if ((($el eq '<')||($el eq '>')||($el eq '-')) && (!$found)){
			$startpos = $counter;
			$found=1;
			}
		$counter++;
		}
	my $left = $startpos-$left_flank;
	my $right = $left+24;
		if ($left <= 0){
		$left_flank--;
		$right_flank++;
		$max_flank--;
		$counter=0;
		next;
		}
		elsif ($right >= $counter){
		$left_flank--;
		$right_flank++;
		$max_flank--;
		$counter=0;
		next;
		}
		else{
		my @structure_line = @struc[$left..$right];
		my @seq_line = @seq[$left..$right];
		my $seqline = join('', @seq_line);
		my $strline = join('', @structure_line);
		push (@structure_lines, $strline);
		push (@sequence_lines, $seqline);
		my $endpos = $startpos+$mir_length-1;
		@mir_struc = @struc[$startpos..$endpos];
		}
	$left_flank--;
	$right_flank++;
	$max_flank--;
	$counter=0;
	}
	my $consecutive_gaps = 0;
	my $failed_gaps = 0;
	my $mirstructure;
		foreach my $basepair (@mir_struc){
		chomp $basepair;
			if ($consecutive_gaps > $max_gaps){
			$failed_gaps = 1;
			}
			if (($basepair ne '.') && ($basepair ne '-')){
			$consecutive_gaps=0;
			}
			else{
			$consecutive_gaps++;
			}
		$mirstructure = "$mirstructure$basepair";
		}
		if ($failed_gaps == 1) {
		$structure = "0";
		my $outseq = 0;
		$mirna_pos = 0;
		return ($structure, $outseq, $mirna_pos);
		}
		else{
		}

	######################################################################################
	# Checks that mirna is involved in base pairing and it is not pairing to itself eg.: #
	# Valid:                                                                             #
	# (((.((((.(((((.......                                                              #
	# ))))))...)))..))))...                                                              #
	# Invalid:                                                                           #
	# .....................                                                              #
	# ((((((....(..))))))))                                                              #
	######################################################################################
	my $basecount = 0;
	
		foreach my $strs (@structure_lines){
		chomp $strs;
		my @no_paired = split(//, $strs);
			foreach my $basepair (@no_paired){
			chomp $basepair;
				if (($basepair ne '.') && ($basepair ne '-')){
				$basecount++;
				}
			}
			#if ((($strs =~m/[\(\<]+/) && ($strs =~m/[\)\>]+/)) || (($strs !~/[\(\<]+/) && ($strs !~/[\)\>]/))){
			if ((($mirstructure =~m/[\(\<]+/) && ($mirstructure =~m/[\)\>]+/)) || (($mirstructure !~/[\(\<]+/) && ($mirstructure !~/[\)\>]/))){
			my $outseq = 0;
			$mirna_pos = 0;
			$structure = "0";
			return ($structure, $outseq, $mirna_pos);
			}
			elsif ($basecount < $min_paired){
			my $outseq = 0;
			$mirna_pos = 0;
			$structure = "0";
			return ($structure, $outseq, $mirna_pos);
			}
		$basecount = 0;
		}
	# This extends the hairpin out as far as it will go..
	if ($bases =~m/[\(\<]/){
	my $done = 0;
	$counter = 0;
		foreach my $el (@mirnaposition){
		chomp $el;
			if (($el=~m/[\<\>\-]/) && (!$done)){
			my $new_count = $counter;
				# Work back to start of hairpin
				while (($struc[$new_count] ne ')') && ($new_count > 0)){
				$new_count--;
				}
				# Could be gaps so work forward until find the first in the hairpin
				while ($struc[$new_count] ne '('){
				$new_count++;
				}
			$seq_start = $new_count;
			$done++;
			}
		$counter++;
		}
		$counter = 0;
		$done = 0;
		my $found =0;
		my $pair_no = $pairs[$seq_start];
		chomp $pair_no;
		# Now need to find where the hairpin ends - we know which bases pair to each so
		# just get the base which pairs with the first base in the hairpin
		foreach my $el (@pairs){
		chomp $el;
			if ($el ne $pair_no){
			$counter++;
			}
			else{
			$seq_end = $counter+1;
			}
		}
	}
	# Same as above but works the other way..
	elsif ($bases =~m/\)/){
	my $done = 0;
	$counter = 0;
		foreach my $el (@mirnaposition){
		chomp $el;
			if (($el ne '<')&&($el ne '>')&&($el ne '-')){
			$counter++;
			}
			elsif (($el=~m/[\<\>\-]/) && (!$done)){
			my $new_count = $counter;
			my $size_count = @struc;
				while (($struc[$new_count] ne '(') && ($new_count < $size_count-1)){
				$new_count++;
				}
				while ($struc[$new_count] ne ')'){
				$new_count--;
				}
			$seq_end = $new_count;
			$done++;
			}
		}
		$counter = 0;
		$done = 0;
		my $found = 0;
		my $pair_no = $pairs[$seq_end];
		chomp $pair_no;
		foreach my $el (@pairs){
		chomp $el;
			if (($el ne $pair_no) && (!$found)){
			$counter++;
			}
			else{
			$seq_start = $counter;
			$found++;
			}
		}
	}
	
my ($outseq, $outstruc, $outmirnapos);
$counter = 0;
# Gets the hairpin sequence
foreach my $el (@seq){
	if (($counter >= $seq_start) && ($counter <= $seq_end)){
	$outseq = "$outseq$el";
	}
	$counter++;
}
$counter =0;
# Gets the hairpin structure
foreach my $el (@struc){
	if (($counter >= $seq_start) && ($counter <= $seq_end)){
	$outstruc = "$outstruc$el";
	}
	$counter++;
}
$structure =$outstruc;
$counter =0;
# Gets the miRNA position in the hairpin
foreach my $el (@mirnaposition){
	if (($counter >= $seq_start) && ($counter <= $seq_end)){
	$outmirnapos = "$outmirnapos$el";
	}
	$counter++;
}
$mirna_pos = $outmirnapos;
chomp $outseq;
my @return;
push (@return, $structure, $outseq, $mirna_pos);
return (@return);
}

# Takes trimmed structure and looks to see if it forms a hairpin and isn't too bulgey

sub analyse_structure{
my $structure = shift;
my $min_hairpin_len = shift ;
my $no_complex_loops = shift;
my $max_unpaired = shift ;
chomp $structure;
my $mirna = 0;
my $struclen = length $structure;
my @structary = split (//,$structure);
	if ($structure =~m/^(\.+)([\(\<]+\.+.+)/){ # Remove leading gaps
	$structure = $2;
	chomp $structure;
	}
	if ($structure =~m/^(.+[\)\>])(\.+)$/){ # Remove trailing gaps
	$structure = $1;
	chomp $structure;
	}
	if (($structure =~m/^[\<\(\.]+\.+[\>\)\.]+$/) && ($struclen >= $min_hairpin_len)){ # Passes
	$mirna = 1;
	}
	elsif (($structure eq 0) || ($struclen < $min_hairpin_len)){ # Fails
	$mirna = 2;
	}
        elsif (($structure =~m/^[\<\(\.]{30,200}.{1,40}[\>\)\.]{30,200}$/) && ($struclen >= $min_hairpin_len)){ # This allows an imperfect hairpin with small sub-structures in the loop region
		if ($no_complex_loops){
		$mirna = 2;
		}
		else{
		$mirna = 3;
		}
        }
	else{
	$mirna = 2;
	}
my $structcount = 0;
my $bulge = 0;
	foreach my $el (@structary){
	$structcount++;
		if (($el eq '.')||($el eq '-')){
		$bulge++;
                }
	}
my $no = ($max_unpaired/100)*$structcount;
	if ($bulge > $no){ # If more than a user defined percentage (or default 33%) of the bases are unpaired we say that it isn't a miRNA
	$mirna = 2;
	}
my $hairpin_length = length($structure);
my @return_analysis;
push (@return_analysis, $mirna);
push (@return_analysis, $hairpin_length);
return (@return_analysis);
}
# Gets GC content of the hairpin
sub get_gc {
  my $dna = shift;
  my $gc = 0;
  my $total = 0;
  for (my $i=0; $i<length($dna); $i=$i+1) {
    if (substr($dna, $i, 1) eq "G" || substr($dna, $i, 1) eq "C") {
      $gc = $gc + 1;
    }
    $total = $total + 1;
  }
  return ($gc/$total)*100;
}

sub find_mir_star{
my $seqstart = shift;
my $seqend = shift;
my $hpseq = shift;
my $hpstruc = shift;
my $mirsequence = shift;
chomp ($hpseq, $hpstruc,$mirsequence);
my $length = $seqend-$seqstart;

chomp ($seqstart, $seqend, $length);
my @hseq = split (//, $hpseq);
my @hstruc = split (//, $hpstruc);
my $i=0;
my $j = @hstruc;
my @pairs = @hstruc;
my $countopen =0;
my @closing;
my @return;

# Find which columns pair with which
	while ($i < $j){
		if (($pairs[$i] eq '(') || ($pairs[$i] eq '<')){
			$countopen++;
			unshift (@closing, $countopen);
			$pairs[$i]=$countopen;
			}
		if (($pairs[$i] eq ')') || ($pairs[$i] eq '>')){
		my $countclose = shift @closing;
		$pairs[$i]=$countclose;
		$countclose--;
		}
	$i++;
	}
$i = 0;
	if ($hpstruc =~m/</){
	my $mirpair;
	my $done = 0;
		foreach my $el (@hstruc){
		chomp $el;
			if (($el eq '<') && (!$done)){
			$mirpair = $pairs[$i];
			chomp $mirpair;
			$done = 1;
			}
			else{
			$i++;
			}
		
		}
	my $found = 0;
	my $pos = -1;
	my $before = 0;
		foreach my $el (@pairs){
		$pos++;
		chomp $el;
			if (($el eq "$mirpair") && (!$found)){
			my $st = $pos;
			my $end = $st+$length;
			my @slice = @hseq[$st..$end];
			my $mir_string;
				foreach my $el (@slice){
				chomp $el;
				$mir_string = "$mir_string$el";
				}
			$found =1;
			}
			elsif (($el eq "$mirpair") && ($found == 1)){
			my $end = $pos+3; # Account for 3' overhang of miRNA*;
				if ($end > (@hseq-1)){
				$end = @hseq-1;
				}
			my $st = ($end-$length)-6; # Account for 3' overhang of miRNA*;
				if ($st < 0){
				$st=0;
				}
			my @slice = @hseq[$st..$end];
			my $star_string;
				foreach my $el (@slice){
				chomp $el;
				$star_string = "$star_string$el";
				}
			my $lf = 0;
				while ($lf < 7){
				my $mirstar = substr($star_string, $lf, $length+1);
				$mirstar =~ tr/U/T/;
				$mirstar =~ tr/a-z/A-Z/;
				push (@return, $mirstar);
				$mirstar = substr($star_string, $lf, $length-1);
				$mirstar =~ tr/U/T/;
				$mirstar =~ tr/a-z/A-Z/;
				push (@return, $mirstar);
				$mirstar = substr($star_string, $lf, $length);
				$mirstar =~ tr/U/T/;
				$mirstar =~ tr/a-z/A-Z/;
				push (@return, $mirstar);
				$lf++;
				}
			return (@return);
			}
		}
	}
	
	else{
	my $done = 0;
	my $mirpair;
	my $p = length($hpstruc)-1;
		for (my $count = $p; $count > 0; $count--) {
		my $el = $hstruc[$count];
		chomp $el;
			if (($el eq '>') && (!$done)){
			$mirpair = $pairs[$count];
			chomp $mirpair;
			$done = 1;
			}
		}
	my $found = 0;
	my $pos = -1;
	my $after = 0;
	my $mpos = 0;
		foreach my $el (@pairs){
		$pos++;
		chomp $el;
			if (($el eq "$mirpair") && (!$found)){
			my $st = $pos-3; # Account for 3' overhang of miRNA*;
				if ($st < 0){
				$st=0;
				}
			my $end = $st+$length+6; # Account for 3' overhang of miRNA*;
				if ($end > @hseq-1){
				$end = @hseq-1;
				}
			my @slice = @hseq[$st..$end];
			my $star_string;
				foreach my $el (@slice){
				chomp $el;
				$star_string = "$star_string$el";
				}
			$found =1;
			my $lf = 0;
				while ($lf < 7){
				my $mirstar = substr($star_string, $lf, $length+1);
				push (@return, $mirstar);
				$lf++;
				}
			}
			elsif (($el eq "$mirpair") && ($found == 1)){
			my $end = $pos;
			my $st = $end-$length;
			my @slice = @hseq[$st..$end];
			my $mir_string;
				foreach my $el (@slice){
				chomp $el;
				$mir_string = "$mir_string$el";
				}
			return (@return);
			}
		}
	}

}
sub validate_mir_star{
my $seqstart = shift;
my $seqend = shift;
my $hpseq = shift;
my $hpstruc = shift;
my $max_gaps = shift;
my $min_paired = shift ;
chomp ($hpseq, $hpstruc);
my $length = $seqend-$seqstart;

chomp ($seqstart, $seqend, $length);
my @hseq = split (//, $hpseq);
my @hstruc = split (//, $hpstruc);
my $i=0;
my $j = @hstruc;
my @pairs = @hstruc;
my $countopen =0;
my @closing;
my $validation = 1;
# Find which columns pair with which
	while ($i < $j){
		if (($pairs[$i] eq '(') || ($pairs[$i] eq '<')){
		$countopen++;
		unshift (@closing, $countopen);
		$pairs[$i]=$countopen;
		}
		elsif (($pairs[$i] eq ')') || ($pairs[$i] eq '>')){
		my $countclose = shift @closing;
		$pairs[$i]=$countclose;
		$countclose--;
		}
	$i++;
	}
$i = 0;
	if ($hpstruc =~m/</){
	my $mirpairstart = 0;
	my $mirpairend = 0;
	my $done = 0;
		foreach my $el (@hstruc){
		chomp $el;
			if (($el eq '<') && (!$done)){
			$mirpairstart = $pairs[$i];
			chomp $mirpairstart;
			$done = 1;
			}
			elsif ($el eq '<'){
			$mirpairend = $pairs[$i];
			chomp $mirpairend;
			}
		$i++;
		}
	my $foundstart = 0;
	my $foundend = 0;
	my $pos = -1;
	my $before = 0;
	my $st = 0;
	my $starstart = 0;
		foreach my $el (@pairs){
		$pos++;
		chomp $el;
			if (($el eq "$mirpairstart") && (!$foundstart)){
			$st = $pos;
			$foundstart = 1;
			}
			elsif (($el eq "$mirpairend") && (!$foundend)){
			my $end = $pos;
			my @slice = @hstruc[$st..$end];
			my @slice2 = @hseq[$st..$end];
			my $mir_string;
			my $mir_seq;
				foreach my $el (@slice){
				chomp $el;
				$mir_string = "$mir_string$el";
				}
				foreach my $el (@slice2){
				chomp $el;
				$mir_seq = "$mir_seq$el";
				}
			$foundend = 1;
			}
			elsif (($el eq "$mirpairend") && ($foundend == 1)){
			$starstart = $pos;
			}
			elsif (($el eq "$mirpairstart") && ($foundstart == 1)){
			my $starend = $pos;
			my @slice = @hstruc[$starstart..$starend];
			my @slice2 = @hseq[$starstart..$starend];
			my $star_string;
			my $star_seq;
				foreach my $el (@slice){
				chomp $el;
				$star_string = "$star_string$el";
				}
				foreach my $el (@slice2){
				chomp $el;
				$star_seq = "$star_seq$el";
				}
			my $lf = 0;
			$validation = &test_mirstar($star_string,$max_gaps,$min_paired);
			return ($validation);
			}
		}
	}
	
	else{
	my $done = 0;
	my $mirpairstart = 0;
	my $mirpairend = 0;
	my $p = length($hpstruc)-1;
		for (my $count = $p; $count > 0; $count--) {
		my $el = $hstruc[$count];
		chomp $el;
			if (($el eq '>') && (!$done)){
			$mirpairend = $pairs[$count];
			chomp $mirpairend;
			$done = 1;
			}
			else{
				if ($el eq '>'){
				$mirpairstart = $pairs[$count];
				chomp $mirpairstart;
				}
			}
		}
	my $foundstart = 0;
	my $foundend = 0;
	my $pos = -1;
	my $after = 0;
	my $mpos = 0;
	my $mirstart = 0;
		foreach my $el (@pairs){
		$pos++;
		chomp $el;
			if (($el eq "$mirpairend") && (!$foundend)){
			$mirstart = $pos;
			$foundend=1;
			}
			elsif (($el eq "$mirpairstart") && (!$foundstart)){
			my $end = $pos;
			my @slice = @hstruc[$mirstart..$end];
			my @slice2 = @hseq[$mirstart..$end];
			my $star_string;
				foreach my $el (@slice){
				chomp $el;
				$star_string = "$star_string$el";
				}
				my $star_seq;
				foreach my $el (@slice2){
				chomp $el;
				$star_seq = "$star_seq$el";
				}
			$foundstart =1;
			my $lf = 0;
			$validation = &test_mirstar($star_string,$max_gaps,$min_paired);
			}
			elsif (($el eq "$mirpairstart") && ($foundstart == 1)){
			$mirstart = $pos;
			}
			elsif (($el eq "$mirpairend") && ($foundend == 1)){
			my $mirend = $pos;
			my @slice = @hstruc[$mirstart..$mirend];
			my @slice2 = @hseq[$mirstart..$mirend];
			my $mir_string;
			my $mir_seq;
				foreach my $el (@slice){
				chomp $el;
				$mir_string = "$mir_string$el";
				}
				foreach my $el (@slice2){
				chomp $el;
				$mir_seq = "$mir_seq$el";
				}
			return ($validation);
			}
		}
	}

}

sub test_mirstar{
my $mirstar = shift;
my $max_gaps = shift ;
my $min_paired = shift ;
my $value = 1;
my $consecutive_unpaired = 0;
my $basecount=0;
my $last = "";
my @no_paired = split(//, $mirstar);
	foreach my $basepair (@no_paired){
	chomp $basepair;
		if ($consecutive_unpaired > $max_gaps){
		$value = 0;
		}
		if (($basepair ne '.') && ($basepair ne '-')){
		$basecount++;
		$consecutive_unpaired=0;
		}
		else{
		$consecutive_unpaired++;
		}
	}
	if ((($mirstar =~m/[\(\<]+/) && ($mirstar =~m/[\)\>]+/)) || (($mirstar !~/[\(\<]+/) && ($mirstar !~/[\)\>]/))){
	$value = 0;
	}
	elsif ($basecount < ($min_paired-2)){
	$value = 0;
	}
return($value);
}
1;