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
# pipeline. This class generates miRNA image files as PDF using
# Ghostscript and ps2pdf..
#
# TODO This class has been converted from a standalone script and might
# still need some cleanup (e.g. passing of parameters)
#
#####################################################################
package SrnaTools::Module::MarkupHairpins ;
use base SrnaTools::Module ;
use strict;
use warnings;
use Cwd;

sub run{
my $self = shift ;

my $cwd = getcwd;
my $workingdir = $self->job->job_working_dir ;
my $input = $workingdir.'/'.$self->param('input_file'); # hairpins.txt

open (INPUT, "$input") or SrnaTools::Exception::ModuleExecution->throw(message=>"failed to open input file",log_msg=>"file: $input");

my ($id, $seq, $structure);
mkdir "$workingdir/data/temp$$" || SrnaTools::Exception::ModuleExecution->throw(message=>"failed to make tmp directory",log_msg=>"tried to make tmp dir in ''$workingdir/data/'");

# parse hairpins.txt
while (<INPUT>){
  if ($_=~m/^>(\S+)\/(\d*\-?\d*)/){ # it's an ID
    $id = $1;
    my $se = $2;
    #print STDERR "$se\n";
    $id = "$id"."_$se";
    chomp $id;
  }
  elsif ($_=~m/^\w+/){ # it's a sequence
    $seq = $_;
    chomp $seq;
    $seq =~ tr/T/U/;
  }
  elsif  ($_=~m/^[\.\(\)\-\<\>\{\}\=]+/){ # it's a structure
    my $nostart = 1;
    my $nostarstart = 1;
    my $miRNAstart = 0;
    my $miRNAend = 0;
    my $miRNAstarstart = 0;
    my $miRNAstarend = 0;
    $structure =$_;
    chomp $structure;
    my @st = split (//, $structure);
    my $count = 0;
    
    # Find start and end pos of miRNA and miRNA*
    foreach my $i (@st){
      $count++;
      if (($nostart) && ($i=~m/[\-\<\>]/)){
        $miRNAstart = $count;
        $nostart =0;
      }
      if ($i=~m/[\-\<\>]/){
        $miRNAend = $count;
      }
      if (($nostarstart) && ($i=~m/[\{\}\=]/)){
        $miRNAstarstart = $count;
        $nostarstart =0;
      }
      if ($i=~m/[\{\}\=]/){
        $miRNAstarend = $count;
      }
    } #foreach
    
    $structure =~ tr/</(/;
    $structure =~ tr/>/)/;
    $structure =~ tr/-/./;
    $structure =~ tr/{/(/;
    $structure =~ tr/}/)/;
    $structure =~ tr/=/./;
    #print STDERR "$seq\n$structure\n";
    open (TMP, ">$workingdir/data/temp$$/$id.tmp") or SrnaTools::Exception::ModuleExecution->throw(message=>"Could not create output file $id.tmp",log_msg=>"full path: $workingdir/data/temp$$/$id.tmp");
    print TMP "$seq\n$structure\n";
    close TMP;
    chdir "$workingdir/data/temp$$/";
    
    if (! -e "$id.tmp"){
      SrnaTools::Exception::ModuleExecution->throw(message=>"File $id.tmp does not exist. Can not run RNAplot",log_msg=>"I am in ".getcwd );
    }
    
    my $rnaplot_cmd ;
    if ($miRNAstarend == 0){
      $rnaplot_cmd="RNAplot --pre \"$miRNAstart $miRNAend 8 0 1 0 omark  1 -10 0\" < $id.tmp" ;
    }
    else{
      $rnaplot_cmd="RNAplot --pre \"$miRNAstart $miRNAend 8 0 1 0 omark $miRNAstarstart $miRNAstarend 8 255 0 255 omark 1 -10 0\" < $id.tmp";
    }
    system($rnaplot_cmd)==0 or SrnaTools::Exception::ModuleExecution->throw(message=>"Could not run RNAplot - administrator needs to check if the binary is available",log_msg=>"tried to run: $rnaplot_cmd");
    
    # rna.ps is automatically created by
    # RNAplot (can't be changed)
    # Copy contents into a new file $id.ps and
    # add a legend 
    open (PS, "rna.ps") or SrnaTools::Exception::ModuleExecution->throw(message=>"Could not open rna.ps");
    
    open (PS_FINAL, ">$workingdir/data/temp$$/$id.ps") or SrnaTools::Exception::ModuleExecution->throw(message=>"Could not open postscript file (ps-final)",log_msg=>"tried to open:$workingdir/data/temp$$/$id.ps "); 
    
    # sneek some legend text into the postscript code
    # NOTE this depends on RNAplot's output of ps
    # so it might break with new version of RNAplot
    while (<PS>) {
      s/%%BeginProlog/\/Helvetica findfont\n12 scalefont\n100 10 moveto\nsetfont\n\($id\) show\n%%BeginProlog/ ;
      print PS_FINAL;
    }
    close PS;
    close PS_FINAL;
    unlink "$id.tmp";
  }
}

chdir "$workingdir/data/temp$$/";
unlink "rna.ps";

# If we have no .ps files now then stop
# This can happen if no candidates were found
# (an empty file hairpins.txt still exists)
my $ps_files = `ls *.ps 2>/dev/null`;
if (!$ps_files){
  return 1;
}

# Run Ghostscript on all files calles *.ps
my $gs_cmd = "gs -sDEVICE=pswrite -sOutputFile=output.ps -dNOPAUSE -dBATCH *.ps 1>/dev/null 2>&1";
if (system($gs_cmd) != 0) {
  SrnaTools::Exception::ModuleExecution->throw(message=>"Could not run Ghostscript - administrator needs to check if the binary is available",log_msg=>"tried to run: $gs_cmd");
}

# ps2pdf creates output.pdf
my $ps2pdf_cmd = "ps2pdf output.ps 1>/dev/null 2>&1";
if (system($ps2pdf_cmd)!=0) {
  SrnaTools::Exception::ModuleExecution->throw(message=>"Could not run ps2pdf - administrator needs to check if the binary is available",log_msg=>"tried to run: $ps2pdf_cmd");
}

my $mv_cmd = "mv output.pdf $workingdir/results/results/structures.pdf";
if (system($mv_cmd)!=0) {
  SrnaTools::Exception::ModuleExecution->throw(message=>"Could not move pdf result file",log_msg=>"tried to run: $mv_cmd");
}
chdir $cwd; # back to where we came from
return 1;
} #run

1;