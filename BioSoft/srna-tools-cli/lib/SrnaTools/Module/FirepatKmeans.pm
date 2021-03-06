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
# It runs the k-means clustering analysis for FiRePat.
# TODO This class has been converted from a standalone script and might
# still need some cleanup 
###

package SrnaTools::Module::FirepatKmeans;
use base SrnaTools::Module ;

#use strict;
#use warnings;

sub run{
my $self = shift ;

my $moduleName = "FiRePat - k-means clustering";

#similarity : Pearson Coefficient
#the permitted number of clusters is the maximum number of k-means calls on the given dataset
#if the permitted number of calls is smaller than the total output of HC the lower part of the suggested series (putative number of clusters) is considered

my $in       = $self->param('corr_infile');
my $clustIn  = $self->param('cluster_infile');
my $output   = $self->param('outfile');
$permittedCls= $self->param('permitted_nr_cls');
$selectCls   = $self->param('selected_nr_cls');

$self->update_status("$moduleName: k-means calculation - this may take a while") ;

open in, $in or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open corr-infile file");

$index = 0;
while(<in>)
{
        chomp;
        @v = split /,/;
        push @data, [@v];
        $timepoints = $#v;
        $index++;
}
close(in);

if($timepoints == 0)
{
	#invalid number of timepoints
	#print "invalid number of timepoints \n";
	$self->job->_run_module(
		'firepat_create_failure_outfiles',
		{
			message => 'Firepat:kmeans:invalid number of timepoints.'
		}
	);
	return 0;
}

if($#data < 0)
{
	#invalid input file
	print "invalid number of timepoints \n";
	$self->job->_run_module(
		'firepat_create_failure_outfiles',
		{
			message => 'Firepat:kmeans:empty input file.'
		}
	);
	return 0;
}

open in, $clustIn or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open cluster infile file");

$indexCluster = 0;
while(<in>)
{       
        chomp;
		if($_ > 0)
		{
			$clsH[$indexCluster] = $_;
			$indexCluster++;
		}
}
close(in);

if($indexCluster == 0)
{
	#invalid number of clusters
	print "invalid number of clusters \n";
	$self->job->_run_module(
		'firepat_create_failure_outfiles',
		{
			message => 'Firepat:kmeans:invalid number of clusters from the hierarchcal clustering method.'
		}
	);
	return 0;
}

if($indexCluster > $permittedCls)
{
        for($i = 0; $i < $permittedCls; $i++)
        {
                $clustersH[$i] = $clsH[$indexCluster - $i];
        }
}
else
{
        for($i = 0; $i < $indexCluster; $i++)
        {
                $clustersH[$i] = $clsH[$i];
        }
        $go = 1;
        while($indexCluster + 1 <= $permittedCls && $go == 1)
        {
                if($indexCluster % 2 == 0)
                {
                        $mid = int(($indexCluster - 1)/ 2)
                }
                else
                {
                        $mid = int(($indexCluster) / 2);
                }
                $new = int(($clustersH[$mid] + $clustersH[$mid+1])/2);
                if($clustersH[$mid] - $new < 2)
                {
                        $go = 0;
                }
                for($i = $indexCluster; $i > $mid; $i--)
                {
                        $clustersH[$i] = $clustersH[$i-1];
                }
                $clustersH[$mid+1] = $new;
                $indexCluster++;
        }
}
if($selectCls == 0)
{
#        $indexCluster--;
}
else
{
        $clustersH[$indexCluster-1] = $selectCls;
}

for($indH = 0; $indH < $indexCluster; $indH++)
{
	#check if there are clusters initialized with 0
#	print "Firepat:KM: number of clusters : $clustersH[$indH] \n";
	if ($clustersH[$indH] != 0)
	{
        $clusters = $clustersH[$indH];
        ###########################################
        $epsilon = $clusters / 1000;

        #choose initial clusters
        for($i = 0; $i < $clusters; $i++)
        {
                $random = int(rand($#data));
                while ($try[$random] != 0)
                {
                        $random = int(rand($#data));
                }
                $try[$random] = 1;
                $init[$i] = $random;
        }

        #initialize centroids
        for($i = 0; $i < $clusters; $i++)
        {
                for($j = 1; $j <= $timepoints; $j++)
                {
                        $centroid[$i][$j] = $data[$init[$i]][$j];
                }
        }

        $change = 1;
        $maxSteps = 100;
        $step = 0;
        while($change != 0 && $step < $maxSteps)
        {
                $change = 0;
                #compute cluster association - distance : Pearson Correlation
                for($i = 0; $i < $index; $i++)
                {
                        $dist[$i]    = 0;
                        $clustId[$i] = -1;
                        for($indexC = 0; $indexC < $clusters; $indexC++)
                        {
                                #compute Pearson
                                $meanX  = 0;
                                $meanY  = 0;
                                $stdevX = 0;
                                $stdevY = 0;
                                $meanXY = 0;
                                for($j = 1; $j <= $timepoints; $j++)
                                {
                                        $meanX  += $data[$i][$j];
                                        $stdevX += $data[$i][$j] * $data[$i][$j];
                                }
                                for($j = 1; $j <= $timepoints; $j++)
                                {
                                        $meanY  += $centroid[$indexC][$j];
                                        $stdevY += $centroid[$indexC][$j] * $centroid[$indexC][$j];
                                }
                                
                                for($j = 1; $j <= $timepoints; $j++)
                                {
                                        $meanXY += $data[$i][$j] * $centroid[$indexC][$j];
                                }
                                
                                $meanX  /= $timepoints;
                                $meanY  /= $timepoints;
                                $meanXY /= $timepoints;
                                $stdevX /= $timepoints;
                                $stdevY /= $timepoints;
                                $stdevX -= $meanX * $meanX;
                                $stdevY -= $meanY * $meanY;
                                $stdevX  = sqrt($stdevX);
                                $stdevY  = sqrt($stdevY);
                                
                                if(($meanY == 0 && $stdevY == 0)||($meanX == 0 && $stdevX == 0))
                                {
                                        $r = 0;
                                }
                                else
                                {
                                        $r = ($meanXY - $meanX * $meanY) / ($stdevX * $stdevY);
                                        
                                        if ($r > $dist[$i])
                                        {
                                                $dist[$i]    = $r;
                                                $clustId[$i] = $indexC;
                                        }
                                }
                        }
                        $cId[$step][$i] = $clustId[$i];
                        
                }#end for
                
                #save old centroids  -  for comparison
                for($i = 0; $i < $clusters; $i++)
                {
                        for($j = 1; $j <= $timepoints; $j++)
                        {
                                $oldCentroid[$i][$j] = $centroid[$i][$j];
                        }
                }
                
                for($i = 0; $i < $clusters; $i++)
                {
                        for($j = 1; $j <= $timepoints; $j++)
                        {
                                $centroid[$i][$j] = 0;
                        }
                }
                
                #compute new centroids
                for($i = 0; $i < $clusters; $i++)
                {
                        $count[$i] = 0;
                }
                
                for($i = 0; $i < $index; $i++)
                {
                        for($j = 1; $j <= $timepoints; $j++)
                        {
                                $centroid[$clustId[$i]][$j] += $data[$i][$j];
                        }
                        $count[$clustId[$i]] ++;
                }
                for($i = 0; $i < $clusters; $i++)
                {
                        if($count[$i] > 0)
                        {
                        for($j = 1; $j <= $timepoints; $j++)
                        {
                                $centroid[$i][$j] /= $count[$i];
                        }
                        }#end if
                }

                for($i = 0; $i < $clusters; $i++)
                {
                        $chg = 0;
                        for($j = 1; $j <= $timepoints; $j++)
                        {
                                $chg += ($centroid[$i][$j] - $oldCentroid[$i][$j]) * ($centroid[$i][$j] - $oldCentroid[$i][$j]);
                        }
                        $chg = sqrt($chg);

                        if($chg  > $epsilon)
                        {
                                $change = 1;
                                last;
                        }
                }
                
                $step++;
        }

        for($i = 0; $i < $index; $i++)
        {
                $memo[$i][0] = $data[$i][0];
                $memo[$i][$indH+1]  = $clustId[$i];
        }
	}#end if clusterH[indH] == 0
}

for($indH = 0; $indH < $indexCluster; $indH++)
{
        $rAvgCls = 0;
        for($ind = 0; $ind < $clustersH[$indH]; $ind++)
        {
        #compute the average intracluster correlation for a given cluster
                $rCls   = 0;
                $addCls = 0;
                for($i = 0; $i <= $index; $i++)
                {
                        for($j = $i+1; $j <= $index; $j++)
                        {
                                if($memo[$i][$indH+1] == $ind && $memo[$j][$indH+1] == $ind)
                                {
                                        #compute correlation
                                        $addCls ++;
                                
                                        $meanX  = 0;
                                        $meanY  = 0;
                                        $stdevX = 0;
                                        $stdevY = 0;
                                        $meanXY = 0;
                                        for($k = 0; $k < $timepoints; $k++)
                                        {
                                                $meanX  += $data[$i][$k];
                                                $stdevX += $data[$i][$k] * $data[$i][$k];
                                        }
                                        for($k = 0; $k < $timepoints; $k++)
                                        {
                                                $meanY  += $data[$j][$k];
                                                $stdevY += $data[$j][$k] * $data[$j][$k];
                                        }
                                        
                                        for($k = 0; $k < $timepoints; $k++)
                                        {
                                                $meanXY += $data[$i][$k] * $data[$j][$k];
                                        }
                                                
                                        $meanX  /= $timepoints;
                                        $meanY  /= $timepoints;
                                        $meanXY /= $timepoints;
                                        $stdevX /= $timepoints;
                                        $stdevY /= $timepoints;
                                        $stdevX -= $meanX * $meanX;
                                        $stdevY -= $meanY * $meanY;
                                        $stdevX  = sqrt($stdevX);
                                        $stdevY  = sqrt($stdevY);
                                        
                                        if(($meanY == 0 && $stdevY == 0)||($meanX == 0 && $stdevX == 0))
                                        {
                                                $r = 0;
                                        }
                                        else
                                        {
                                                $r = ($meanXY - $meanX * $meanY) / ($stdevX * $stdevY);
                                        }
                                        
                                        $rCls += $r;
                                }
                        }
                }
                if($addCls == 0)
                {
                        $rCls = 0;
                }
                else
                {
                        $rCls /= $addCls;
                }
        
                $rAvgCls += $rCls;
        }
        
        $rAvgCls /= $clustersH[$indH];
        $intra[$indH] = $rAvgCls;
        
        for($i = 0; $i < $clustersH[$indH]; $i++)
        {
                for($j = 0; $j < $timepoints; $j++)
                {
                        $centroid[$i][$j] = 0;
                }
                $add[$i] = 0;
        }
        for($i = 0; $i < $index; $i++)
        {
                $add[$memo[$i][$indH+1]]++;
                for($j = 0; $j < $timepoints; $j++)
                {
                        $centroid[$memo[$i][$indH+1]][$j] += $data[$i][$j+1];
                }
        }
        for($i = 0; $i < $clustersH[$indH]; $i++)
        {
#               if($add[$i] == 0)
#               {
#                       print "$i : add : $add[$i] \n";
#               }
                for($j = 0; $j < $timepoints; $j++)
                {
                        if($add[$i] == 0)
                        {
                                $centroid[$i][$j] = 0;
                        }
                        else
                        {
                                $centroid[$i][$j] /= $add[$i];
                        }
                }
        }
        
        $rInterAvg = 0;
        for($i = 0; $i < $clustersH[$indH]; $i++)
        {
                for($j = $i+1; $j < $clustersH[$indH]; $j++)
                {
                        #compute Pearson
                        $meanX  = 0;
                        $meanY  = 0;
                        $stdevX = 0;
                        $stdevY = 0;
                        $meanXY = 0;
                        for($k = 0; $k < $timepoints; $k++)
                        {
                                $meanX  += $centroid[$i][$k];
                                $stdevX += $centroid[$i][$k] * $centroid[$i][$k];
                        }
                        for($k = 0; $k < $timepoints; $k++)
                        {
                                $meanY  += $centroid[$j][$k];
                                $stdevY += $centroid[$j][$k] * $centroid[$j][$k];
                        }
                        
                        for($k = 0; $k < $timepoints; $k++)
                        {
                                $meanXY += $centroid[$i][$k] * $centroid[$j][$k];
                        }
                                
                        $meanX  /= $timepoints;
                        $meanY  /= $timepoints;
                        $meanXY /= $timepoints;
                        $stdevX /= $timepoints;
                        $stdevY /= $timepoints;
                        $stdevX -= $meanX * $meanX;
                        $stdevY -= $meanY * $meanY;
                        $stdevX  = sqrt($stdevX);
                        $stdevY  = sqrt($stdevY);
                        
                        if(($meanY == 0 && $stdevY == 0)||($meanX == 0 && $stdevX == 0))
                        {
                                $r = 0;
                        }
                        else
                        {
                                $r = ($meanXY - $meanX * $meanY) / ($stdevX * $stdevY);
                        }
                        $rInterAvg += $r;
                }
        }
        if($clustersH[$indH] > 1)
        {
                $rInterAvg /= $clustersH[$indH] * ($clustersH[$indH] - 1) / 2;
        }
        $inter[$indH] = $rInterAvg;
}

$maxPoz = -1;$maxScore = 0;
for($indH = 0; $indH < $indexCluster; $indH++)
{
        $score = 0.4 * $intra[$indH] + 0.6 * (1- abs($inter[$indH]));
        if($score > $maxScore)
        {
                $maxScore = $score;
                $maxPoz = $indH;
        }
}

close(out);
open out, ">",$output or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output file");

if($selectCls != 0)
{
        $maxPoz = $indexCluster-1;
}       

for($i = 0; $i < $index; $i++)
{
        print out "$memo[$i][0], $memo[$i][$maxPoz + 1] \n";
		#print "$memo[$i][0], $memo[$i][$maxPoz + 1] \n";
}

return 1;
}

## Declare the subroutines
sub trim($);

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

return 1;