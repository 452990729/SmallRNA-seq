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
# TODO This class has been converted from a standalone script and might
# still need some cleanup 
###

package SrnaTools::Module::FirepatCreateHtml;
use base SrnaTools::Module ;

#use strict;
#use warnings;

sub run{
my $self = shift ;

my $moduleName = "FiRePat - create html";
##create html document
#acording to options it can show
# 1. sRNA time series
# 2. microarray time series
# 3. BLAST informations
# 4. similarity information (Pearson Coefficient)
# 5. clusters information
# and combinations of all of the above, as described in the config file
# use :           perl create_html.pl 12345 config intervals
$options    = $self->param('options');
$intervals  = $self->param('color_int');
$fileLoc    = $self->param('file_loc');
$output     = $self->param('res_loc');
#$output     = $self->job->job_working_dir ;

$self->update_status("$moduleName: creating final html") ;

my @config = (
'sRNA_time_series corr_in1L.csv corr_in1C.csv',
'microarray_time_series corr_in2L.csv corr_in2C.csv',
'similarity_information corr.csv',
'clusters_information OUTcls.csv',
);

$options .= '#';
@opt = split(//, $options);

#read the config file
$i = 0;

#flags for sorting the output
$similarity = 0;
$clusters   = 0;

foreach (@config)
{
        chomp;
        @option = split (/ /,$_);
        my $case=$opt[$i];
        { 
                if ($case==1)
                {
					#print "loci data : $fileLoc$option[1] \n";
                        open Lnum, $fileLoc.$option[1] or SrnaTools::Exception::ModuleExecution->throw( "cannot open loci numeric file");
                        open Lcol, $fileLoc.$option[2] or SrnaTools::Exception::ModuleExecution->throw( "cannot open loci color file");
                        $#lociNum = -1;
						$#lociCol = -1;
                        #read numerical values and colors
                        while(<Lnum>)
                        {
                                chomp;
                                @num = split /,/;
                                $timepoints = $#num - 1;
                                push @lociNum, [@num];
                        }
						close(Lnum);
                        while(<Lcol>)
                        {
                                chomp;
                                @col = split /,/;
                                push @lociCol, [@col];
                        }
						close(Lcol);
						#print "loci numeric : $#lociNum \n";
						#print "loci color : $#lociCol \n";
						if($#lociNum < 0 || $#lociCol < 0)
						{
							$self->job->_run_module(
								'firepat_create_failure_outfiles',
								{
									message => 'Firepat:createHTML:empty input file 1.'
								}
							);
							return 0;
						}
                }
                elsif ($case==2)
                {
						#print "microarray : $fileLoc$option[1]\n";
                        open Mnum, $fileLoc.$option[1] or SrnaTools::Exception::ModuleExecution->throw( "cannot open microarray numeric file");
                        open Mcol, $fileLoc.$option[2] or SrnaTools::Exception::ModuleExecution->throw( "cannot open microarray color file");
                        $#microNum = -1;
						$#microCol = -1;
                        #read numerical values and colors
                        while(<Mnum>)
                        {
                                chomp;
								
                                @num = split /,/;
                                $timepoints = $#num - 1;
                                #$timepoints = $#num;
                                push @microNum, [@num];
                        }
						close(Mnum);
                        while(<Mcol>)
                        {
                                chomp;
								
                                @col = split /,/;
                                push @microCol, [@col];
                        }
						close(Mcol);
						#print "microarray numberic: $#microNum \n";
						#print "microarray color : $#microCol \n";
						if($#microNum < 0 || $#microCol < 0)
						{
							$self->job->_run_module(
								'firepat_create_failure_outfiles',
								{
									message => 'Firepat:createHTML:empty input file 2.'
								}
							);
							return 0;
						}
                }
                elsif ($case==3)
                {
                } 
				elsif ($case==4)
                {
                        $similarity = 1;
                        open sim, $fileLoc.$option[1] or SrnaTools::Exception::ModuleExecution->throw( "cannot open similarity file");
                        $#similar = -1;
						while(<sim>)
                        {
                            chomp;
							
                            @s = split /,/;
                            push @similar, [(@s)];
                        }
						close(sim);
						#print "similairty pairs : $#similar \n";
						
                        if($#similar < 0)
						{
							$self->job->_run_module(
								'firepat_create_failure_outfiles',
								{
									message => 'Firepat:createHTML:empty similarity input.'
								}
							);
							return 0;
						}
                }
                elsif ($case==5)
                {
					#print "check clusters \n";
                        $clusters = 1;
						
						$fileName = $fileLoc.$option[1];
						
                        open clsOnLoci, $fileName or SrnaTools::Exception::ModuleExecution->throw( "cannot open input cluster file $fileName");
						system("ls -l $fileLoc > /dev/null");
						
						$#cluster = -1;
						chomp;
                        while(<clsOnLoci>)
                        {
                                chomp;
                                @clstr = split/,/;
                                push @cluster, [@clstr];
                        }
						close(clsOnLoci);
						
						#print "cls count : $#cluster \n";
						if($#cluster < 0)
						{
							$self->job->_run_module(
								'firepat_create_failure_outfiles',
								{
									message => 'Firepat:createHTML:empty cluster input.'
								}
							);
							return 0;
						}
                }
        }
        $i++;
}
if($#microNum > $#lociNum)
{
        $dataLen = $#microNum;
}
else
{
        $dataLen = $#lociNum;
}

#initialize the order
for($i = 0; $i < $dataLen; $i++)
{
        $order[$i] = $i;
}

if ($clusters   == 1)
{
        #sort by $cluster[$i][1]
        $sorted = 0;
        $index = 0;
        while ($sorted == 0 && $index < 100000)
        {
                $sorted = 1;
                for($i = 0; $i < $dataLen - 1; $i++)
                {
                        if($cluster[$order[$i]][1] > $cluster[$order[$i + 1]][1])
                        {
                                $aux         = $order[$i];
                                $order[$i]   = $order[$i+1];
                                $order[$i+1] = $aux;
                                $sorted = 0;
                        }#end if
                }#end for
                $index++;
        }#end while
}#end if

#open outP, '>',$output.'Poz.html' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output html file");
open outP, '>',$output.'/positively_correlated.html' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output html file");

#open outN, '>',$output.'Neg.html' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output html (neg) file");
open outN, '>',$output.'/negatively_correlated.html' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output html (neg) file");

$nameFile = $fileLoc.'names.txt';
open names, $nameFile or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open names file");

while(<names>)
{
        @n = split/,/;
        push @names,[@n];
}

close names;
        
print outP "<html>\n";
print outP "<head>\n";
print outP "<title>Microarray & sRNA analysis </title>\n";
print outP "</head>\n";
print outP "<body>\n";

print outP "<table> \n";

print outN "<html>\n";
print outN "<head>\n";
print outN "<title>Microarray & sRNA analysis </title>\n";
print outN "</head>\n";
print outN "<body>\n";

print outN "<table> \n";

$intLength = 200 / $intervals;
$prev = -1;

print outP "number of timepoints : $#n \n";
print outP "selected pairs : $dataLen \n";
print outP "\n";
print outP "<tr>";
for($i = 0; $i <= $#n; $i++)
{
        print outP "<td>$names[0][$i]</td>";
}
for($i = 0; $i <= $#n; $i++)
{
        print outP "<td>$names[1][$i]</td>";
}
print  outP "<td>Correlation</td>";
print  outP "<td>ClusterId</td><tr>\n";


print outN "number of timepoints : $#n \n";
print outN "selected pairs : $dataLen \n";
print outN "\n";
print outN "<tr>";
for($i = 0; $i <= $#n; $i++)
{
        print outN "<td>$names[0][$i]</td>";
}
for($i = 0; $i <= $#n; $i++)
{
        print outN "<td>$names[1][$i]</td>";
}
print  outN "<td>Correlation</td>";
print  outN "<td>ClusterId</td><tr>\n";

#print to colored html file
for($i = 0; $i < $dataLen; $i++)
{
        
if ($i > 0 && $clusters == 1 && $cluster[$order[$i]][1] != $prev)
{
        if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
        {
                print outP "<tr><td></td> <td></td> <td></td> <td></td> </tr> \n";
                print outP "<tr><td></td> <td></td> <td></td> <td></td> </tr> \n";
        }
        else
        {
                print outN "<tr><td></td> <td></td> <td></td> <td></td> </tr> \n";
                print outN "<tr><td></td> <td></td> <td></td> <td></td> </tr> \n";
        }
}
        if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
        {
                print outP "<tr> \n";
        }
        else
        {
                print outN "<tr> \n";
        }
        for($k = 0; $k < $#opt; $k++)
        {
                my $case = $opt[$k];
                {
                        if ($case==1)
                        {
                                @name  = split(/\//,$lociNum[$order[$i]][0]);
                                @limit = split(/-/,$name[1]);
                                #################################################################################
                                #  change the link format
                                #################################################################################
                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                {
                                        print outP "<td>$lociNum[$order[$i]][0]</td> \n";
                                }
                                else
                                {
                                        print outN "<td>$lociNum[$order[$i]][0]</td> \n";
                                }
                                for($j = 1; $j <= $timepoints; $j++)
                                {
                                        if ($lociCol[$order[$i]][2*$j - 1] == 0 && $lociCol[$order[$i]][2*$j] == 0)
                                        {
                                                $red   = 0;
                                                $green = 0;
                                        }
                                        else
                                        {
                                                $red   = 0;
                                                $green = 0;
                                                if($lociCol[$order[$i]][2*$j - 1] != 0)
                                                {
                                                        $red   = $intLength * $lociCol[$order[$i]][2*$j - 1] + 50;
                                                }
                                                if ($lociCol[$order[$i]][2*$j] != 0)
                                                {
                                                        $green = $intLength * $lociCol[$order[$i]][2*$j]     + 50;
                                                }
                                        }
                                        
                                        if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                        {
                                                printf outP "<td style=\"background-color:rgb($red, $green,0)\" BGCOLOR=\"rgb($red, $green,0)\">";
                                        }
                                        else
                                        {
                                                printf outN "<td style=\"background-color:rgb($red, $green,0)\" BGCOLOR=\"rgb($red, $green,0)\">";
                                        }
                                        
                                        if($red < 125 || $green < 125)
                                        {
                                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                                {
                                                        print outP "<font color=\"white\">";
                                                        printf outP "%0.3f ",$lociNum[$order[$i]][$j];
                                                        print outP "</font>";
                                                }
                                                else
                                                {
                                                        print outN "<font color=\"white\">";
                                                        printf outN "%0.3f ",$lociNum[$order[$i]][$j];
                                                        print outN "</font>";
                                                }
                                        }
                                        else
                                        {
                                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                                {
                                                        print outP "%0.3f",$lociNum[$order[$i]][$j];
                                                }
                                                else
                                                {
                                                        print outN "%0.3f",$lociNum[$order[$i]][$j];
                                                }
                                        }
                                        if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                        {
                                                print outP "</td> \n";
                                        }
                                        else
                                        {
                                                print outN "</td> \n";
                                        }
                                }#end for j
                        }
                        elsif ($case == 2)
                        {
                        
                                #################################################################################
                                #  change the link format -- insert link for arabidopsis NO!
                                #################################################################################
                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                {
                                        print outP "<td>$microNum[$order[$i]][0]</td> \n";
                                }
                                else
                                {
                                        print outN "<td>$microNum[$order[$i]][0]</td> \n";
                                }
                                for($j = 1; $j <= $timepoints; $j++)
                                {
                                        if ($microCol[$order[$i]][2*$j - 1] == 0 && $microCol[$order[$i]][2*$j] == 0)
                                        {
                                                $red   = 0;
                                                $green = 0;
                                        }
                                        else
                                        {
                                                $red   = 0;
                                                $green = 0;
                                                if($microCol[$order[$i]][2*$j - 1] != 0)
                                                {
                                                        $red   = $intLength * $microCol[$order[$i]][2*$j - 1] + 50;
                                                }
                                                if ($microCol[$order[$i]][2*$j] != 0)
                                                {
                                                        $green = $intLength * $microCol[$order[$i]][2*$j]     + 50;
                                                }
                                        }
                                        if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                        {
                                                printf outP "<td style=\"background-color:rgb($red, $green,0)\" BGCOLOR=\"rgb($red, $green,0)\">";
                                                if ($red < 125 || $green < 125)
                                                {
                                                        print outP "<font color=\"white\">";
                                                        printf outP "%0.3f ",$microNum[$order[$i]][$j];
                                                        print outP "</font>";
                                                }
                                                else
                                                {
                                                        printf outP "%0.3f ",$microNum[$order[$i]][$j];
                                                }
                                                printf outP "</td> \n";
                                        }
                                        else
                                        {
                                                printf outN "<td style=\"background-color:rgb($red, $green,0)\" BGCOLOR=\"rgb($red, $green,0)\">";
                                                if ($red < 125 || $green < 125)
                                                {
                                                        print outN "<font color=\"white\">";
                                                        printf outN "%0.3f ",$microNum[$order[$i]][$j];
                                                        print outN "</font>";
                                                }
                                                else
                                                {
                                                        printf outN "%0.3f ",$microNum[$order[$i]][$j];
                                                }
                                                printf outN "</td> \n";
                                        }
                                }

                        }
                        elsif ($case == 3)
                        {
                                #blast information
                        }
                        elsif ($case == 4)
                        {
                                #similarity information
                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                {
                                        print  outP "<td>";
                                        printf outP "%0.3f ", $similar[$order[$i]][2 * ($timepoints + 1)];
                                        print  outP "</td>\n";
                                }
                                else
                                {
                                        print  outN "<td>";
                                        printf outN "%0.3f ", $similar[$order[$i]][2 * ($timepoints + 1)];
                                        print  outN "</td>\n";

                                }
                        }
                        elsif ($case == 5)
                        {
                                #cluster information
                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                {
                                        print  outP "<td>";
                                        printf outP " $cluster[$order[$i]][1]";
                                        print  outP "</td>\n";
                                }
                                else
                                {
                                        print  outN "<td>";
                                        printf outN " $cluster[$order[$i]][1]";
                                        print  outN "</td>\n";

                                }
                        }
                }
        }#end for k - options
        if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
        {
                print outP "</tr>\n";
        }
        else
        {
                print outN "</tr>\n";
        }

        $prev =  $cluster[$order[$i]][1];
#}#end if similarity ... check for better version!!
#}#end if cluster != -1
}#end for i
print outP "</table>\n";
print outP "</body>\n";
print outP "</html>\n";
close(outP);

print outN "</table>\n";
print outN "</body>\n";
print outN "</html>\n";
close(outN);


#=================================================

open outP, '>',$output.'/positively_correlated.csv' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output csv file");
open outN, '>',$output.'/negatively_correlated.csv' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output csv (neg) file");

#open outP, '>',$output.'Poz.csv' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output csv file");
#open outN, '>',$output.'Neg.csv' or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output csv (neg) file");

for($i = 0; $i <= $#n; $i++)
{
        print outP "$names[0][$i],";
}
for($i = 0; $i <= $#n; $i++)
{
        print outP "$names[1][$i],";
}
print  outP "Correlation,";
print  outP "ClusterId\n";

for($i = 0; $i <= $#n; $i++)
{
        print outN "$names[0][$i],";
}
for($i = 0; $i <= $#n; $i++)
{
        print outN "$names[1][$i],";
}
print  outN "Correlation,";
print  outN "ClusterId\n";


#print to csv file
for($i = 0; $i < $dataLen; $i++)
{
        
        for($k = 0; $k < $#opt; $k++)
        {
                my $case = $opt[$k] ;
                {
                        if ($case == 1)
                        {
                                @name  = split(/\//,$lociNum[$order[$i]][0]);
                                @limit = split(/-/,$name[1]);
                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                {
                                        print outP "$lociNum[$order[$i]][0],";
                                }
                                else
                                {
                                        print outN "$lociNum[$order[$i]][0],";
                                }
                                for($j = 1; $j <= $timepoints; $j++)
                                {
                                        if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                        {
                                                printf outP "%0.3f, ",$lociNum[$order[$i]][$j];
                                        }
                                        else
                                        {
                                                printf outN "%0.3f, ",$lociNum[$order[$i]][$j];
                                        }
                                }#end for j
                        }
                        elsif ($case == 2)
                        {
                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                {
                                        print outP "$microNum[$order[$i]][0],";
                                }
                                else
                                {
                                        print outN "$microNum[$order[$i]][0],";
                                }
                                
                                
                                for($j = 1; $j <= $timepoints; $j++)
                                {
                                        if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                        {
                                                printf outP "%0.3f, ",$microNum[$order[$i]][$j];
                                        }
                                        else
                                        {
                                                printf outN "%0.3f, ",$microNum[$order[$i]][$j];
                                        }
                                }
                        }
                        elsif ($case == 3)
                        {
                                #blast information
                        }
                        elsif ($case == 4)
                        {
                                #similarity information
                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                {
                                        printf outP "%0.3f, ", $similar[$order[$i]][2 * ($timepoints + 1)];
                                }
                                else
                                {
                                        printf outN "%0.3f, ", $similar[$order[$i]][2 * ($timepoints + 1)];
                                }
                        }
                        elsif ($case == 5)
                        {
                                #cluster information
                                if($similar[$order[$i]][2 * ($timepoints + 1)] > 0)
                                {
                                        printf outP " $cluster[$order[$i]][1] \n";
                                }
                                else
                                {
                                        printf outN " $cluster[$order[$i]][1] \n";
                                }
                        }
                }
        }#end for k - options

        $prev =  $cluster[$order[$i]][1];
#}#end if similarity ... check for better version!!
#}#end if cluster != -1
}#end for i
return 1;
}

return 1;