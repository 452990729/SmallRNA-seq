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
# It calculates log2 ratios for FiRePat.
# TODO This class has been converted from a standalone script and might
# still need some cleanup 
###

package SrnaTools::Module::FirepatCreateColors;
use base SrnaTools::Module ;

#use strict;
#use warnings;

sub run{
my $self = shift ;

my $moduleName = "FiRePat - coloring";
#compute log2 ratios relative to first timepoint
$input = $self->param('infile');
$output= $self->param('outfile');
$intervals = $self->param('color_int');

open in, $input or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open input file file");
open out, '>', $output or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output file file");

# creates the color pattern for the input file
# the input file should be csv file
#    perl create_colors.pl  input_numeric output_color intervals
# uses normal distribution of colors.
# green = pozitive values
# red   - negative values

$#micro = -1;
#$Zthreshold = 3;
$Zthreshold = 1;
$index = 0;
while(<in>)
{
        chomp;
        @data = split /,/;
        $timepoints = $#data - 1;
        push @micro, [@data];
        $index++;
}#end while
close(in);

if($#micro < 0)
{
	$self->job->_run_module(
		'firepat_create_failure_outfiles',
		{
			message => 'Firepat:CreateColors:empty input file.'
		}
	);
	return 0;
}

#adjust Z threshold for the number of entries
if($index < 100)
{
        $Zthreshold = 1;
}
else
{
        if($index < 200)
        {
                $Zthreshold = 2;
        }
        else
        {
                $Zthreshold = 3;
        }
}

#print "index = $index \n";
#print "timepoints = $timepoints \n";
#compute mean and stdev
$mean  = 0;
$stdev = 0;
$count = 0;
for($i = 0; $i < $index; $i++)
{
        for($j = 1; $j <= $timepoints; $j++)
        {
                if($micro[$i][$j] != 0)
                {
                        $mean  += $micro[$i][$j];
                        $stdev += $micro[$i][$j] * $micro[$i][$j];
                        $count++;
                }
        }
}
if($count > 0)
{
	$mean  /= $count;
	$stdev /= $count;
	$stdev -= $mean * $mean;
	$stdev  = sqrt($stdev);
}
else
{
	$mean = 0;
	$stdev = 0;
}

#compute Z scores
for($i = 0; $i < $index; $i++)
{
        for($j = 0; $j < $timepoints; $j++)
        {
			if($stdev > 0)
			{
                $z[$i][$j] = ($micro[$i][$j] - $mean) / $stdev;
			}
			else
			{
				$z[$i][$j] = 0;
			}
        }
}

#compute intervals
#compute the limits
$epsU = 1;
$epsD = 1;
$memoIU = -1;
$memoJU = -1;
$memoID = -1;
$memoJD = -1;

for($i = 0; $i < $index; $i++)
{
        for($j = 1; $j <= $timepoints; $j++)
        {
                if ($z[$i][$j] >= 0)
                {
                        if(abs($z[$i][$j] - $Zthreshold) < $epsU)
                        {
                                $zMax = $z[$i][$j];
                                $memoIU = $i;
                                $memoJU = $j;
                                $epsU = abs($z[$i][$j] - $Zthreshold);
                        }
                }
                else
                {
                        #print "$z[$i][$j] ";
                        if(abs($z[$i][$j] + $Zthreshold) < $epsD)
                        {
                                $zMin = $z[$i][$j];
                                $memoID = $i;
                                $memoJD = $j;
                                $epsD = abs($z[$i][$j] + $Zthreshold);
                        }
                }
        }
}

# print "zMax = $zMax \n";
# print "microMax = $micro[$memoIU][$memoJU]\n";
# print "zMin = $zMin \n";
# print "microMin = $micro[$memoID][$memoJD]\n";
# print "number of intervals [on the pozitive side and on the negative side] : $intervals \n";
#interval length
$iLenP = 0.65 * $micro[$memoIU][$memoJU] / $intervals;
$iLenN = 0.65 * (- $micro[$memoID][$memoJD]) / $intervals;
# print "interval length pozitive: $iLenP \n";
# print "interval length negative: $iLenN \n";
#equally divided intervals between min and 0  and 0 and max
for($i = 0; $i <= $intervals; $i++)
{
        $intP[$i] = $i * $iLenP;
}

for($i = 0; $i <= $intervals; $i++)
{
        $intN[$i] = $i * $iLenN + 0.65 * $micro[$memoID][$memoJD] ;
}
#associate colors
for($i = 0; $i < $index; $i++)
{
        for($j = 1; $j <= $timepoints; $j++)
        {
                $green = 0;
                $red   = 0;
                if ($micro[$i][$j] >= 0)
                {
                        #green
                        $interval = 0;
                        while ($intP[$interval] < $micro[$i][$j] && $interval <= $intervals)
                        {
                                $interval++;
                        }
                        $green = $interval;
                }
                else
                {
                        #red
                        $interval = 0;
                        while ($intN[$interval] < $micro[$i][$j] && $interval <= $intervals)
                        {
                                $interval++;
                        }
                        $red = $intervals - $interval;
                }
                $color[$i][$j][0] = $red;
                $color[$i][$j][1] = $green;
        }#end for j
}#end for i

#print color data to file
#for each object it prints the id
#and for each timepoint the red and the green associated interval
for($i = 0; $i < $index; $i++)
{
        print out "$micro[$i][0],";
        for($j = 1; $j <= $timepoints; $j++)
        {
                print out "$color[$i][$j][0], $color[$i][$j][1],";
        }
        print out "\n";
}
close(out);

return 1;
}

return 1;