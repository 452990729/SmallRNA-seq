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
# It runs the clustering analysis for FiRePat.
# TODO This class has been converted from a standalone script and might
# still need some cleanup (e.g. passing of parameters)
###

package SrnaTools::Module::FirepatHierarchicalCluster;
use base SrnaTools::Module ;

#use strict;
#use warnings;

sub run{
my $self = shift ;

$in = '';
$clusterThreshold = -1;
$outfile = '';

my $moduleName = "FiRePat - hierarchical clustering";
#input sRNA file : definition of the loci
$in = $self->param('infile');
$clusterThreshold = $self->param('nr_cluster_threshold');
$outfile = $self->param('outfile');

$self->update_status("$moduleName: finding clusters - this may take a while") ;

#maximum number of levels for the computed dendrogram
$maxStep = 100;
open in, $in or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open input file");

$#data = -1;
while(<in>)
{
        chomp;
        @v = split /,/;
        $timepoints = $#v;
        push @data, [@v];
        $index++;
}
close(in);

#check if file empty
if($#data < 0 || $timepoints == 0)
{
#	print "Firepat:HC:no timepoints or no data in the input file \n";
	#input file does not exist
	$self->job->_run_module(
		'firepat_create_failure_outfiles',
		{
			message => 'empty input file.'
		}
	);
	return 0;
}
#print "Firepat:HC: input file ok \n";

$step = 0;
#initialize the clusters
for($i = 0; $i < $#data; $i++)
{
#  $self->update_status("$moduleName: finding clusters - this may take a while (initialising step $i of ".$#data.")") ;
  
        for($j = 1; $j <= $timepoints; $j++)
        {
                $centroids[$step][$i][$j - 1] = $data[$i][$j];
        }
}
$clusters[$step] = $#data;

#each line will contain the number of clusters
open out, '>',$outfile or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output file");

while ($clusters[$step] > $clusterThreshold && $step < $maxStep)
{

#  $self->update_status("$moduleName: finding clusters - this may take a while (clustering step $step of a maximum ".$maxStep.")") ;
  
        # compute Pearson correlation and create the association table
        # (most_similar_cluster, similarity)
        for($i = 0; $i < $clusters[$step]; $i++)
        {
                $cluster[$i][0] = -1; # most similar neighbour
                $cluster[$i][1] =  $pearsonThreshold; # closest distance
                for($j = 0; $j < $clusters[$step]; $j++)
                {
                        $meanX  = 0;
                        $meanY  = 0;
                        $meanXY = 0;
                        $stdevX = 0;
                        $stdevY = 0;
                        $r      = 0;
                        for($k = 1; $k <= $timepoints; $k++)
                        {
                                $meanX  += $centroids[$step][$i][$k];
                                $stdevX += $centroids[$step][$i][$k] * $centroids[$step][$i][$k];
                                
                                $meanY  += $centroids[$step][$j][$k];
                                $stdevY += $centroids[$step][$j][$k] * $centroids[$step][$j][$k];
                                
                                $meanXY += $centroids[$step][$i][$k] * $centroids[$step][$j][$k];
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
                        
                        if($stdevX == 0 || $stdevY == 0)
                        {
                                $r = 0;
                        }
                        else
                        {
                                $r = ($meanXY - $meanX * $meanY) / ($stdevX * $stdevY);
                        }
                        if ($r > $cluster[$i][1] && $i != $j)
                        {
                                $cluster[$i][0] = $j;
                                $cluster[$i][1] = $r;
                        }
                }
        }
        
        #assign clusters
        $nrCls = 1;
        #initialize clusters
        for($i = 0; $i < $clusters[$step]; $i++)
        {
                $cls[$i] = 0;
        }
        
        for($i = 0; $i < $clusters[$step]; $i++)
        {
                if($cluster[$i][0] == -1)
                {
                        $cls[$i] = $nrCls;
                        $nrCls++;
                        next;
                }
                if($cls[$i] == 0 && $cls[$cluster[$i][0]] == 0)
                {
                        $cls[$i] = $nrCls;
                        $cls[$cluster[$i][0]] = $nrCls;
                        $nrCls++;
                        next;
                }
                if($cls[$i] == 0 && $cls[$cluster[$i][0]] != 0)
                {
                        $cls[$i] = $cls[$cluster[$i][0]];
                        next;
                }
                if($cls[$i] != 0 && $cls[$cluster[$i][0]] == 0)
                {
                        $cls[$cluster[$i][0]] = $cls[$i];
                        $next;
                }
                if($cls[$i] != 0 && $cls[$cluster[$i][0]] != 0)
                {
                        if ($cls[$i] < $cls[$cluster[$i][0]])
                        {
                                for($j = 0; $j <=$clusters[$step]; $j++)
                                {
                                        if ($cls[$j] == $cls[$cluster[$i][0]])
                                        {
                                                $cls[$j] = $cls[$i];
                                        }
                                }
                                $cls[$cluster[$i][0]] = $cls[$i];
                        }
                        else
                        {
                                for($j = 0; $j <=$clusters[$step]; $j++)
                                {
                                        if ($cls[$j] == $cls[$i])
                                        {
                                                $cls[$j] = $cls[$cluster[$i][0]];
                                        }
                                }
                                $cls[$i] = $cls[$cluster[$i][0]];
                        }       
                }
        }
        
        for($i = 0; $i < $clusters[$step]; $i++)
        {
                $missing[$i] = 0;
        }
        for($i = 0; $i < $clusters[$step]; $i++)
        {
                $missing[$cls[$i]] = 1;
        }
        
        $nrClsNew = $nrCls;
        for($i = $nrCls - 1; $i > 0; $i--)
        {
                if($missing[$i] == 0)
                {
                        $nrClsNew --;
                        for($j = 0; $j < $clusters[$step]; $j++)
                        {
                                if ($cls[$j] > $i)
                                {
                                        $cls[$j]--;
                                }
                        }
                }
        }

        $nrClsNew--;
        
        print out "$nrClsNew\n";

        for($i = 0; $i <= $nrClsNew; $i++)
        {
                $count[$i] = 0;
        }
        #compute centroids
        for($i = 0; $i < $nrClsNew; $i++)
        {
                for($j = 0; $j < $timepoints; $j++)
                {
                        $cen[$i][$j] = 0;
                }
        }
        for($i = 0; $i < $clusters[$step]; $i++)
        {
                for($j = 0; $j < $timepoints; $j++)
                {
                        $cen[$cls[$i]][$j] += $centroids[$step][$i][$j]; 
                }
                $count[$cls[$i]]++;
        }
        
        for($i = 1; $i <= $nrClsNew; $i++)
        {
                for($j = 0; $j < $timepoints; $j++)
                {
					if ($count[$i] > 0)
					{
                        $cen[$i][$j] /= $count[$i];
					}
					else
					{
						$cen[$i][$j] = 0;
					}
                }
        }
        
        $step++;
        #shift centroids : clusters 1 - <= nrCls => 0 .. < nrCls 
        $clusters[$step] = $nrClsNew;
        for($i = 0; $i < $clusters[$step]; $i++)
        {
                for($j = 0; $j < $timepoints; $j++)
                {
                        $centroids[$step][$i][$j] = $cen[$i + 1][$j];
                }
        }
        
        if($clusters[$step] == $clusters[$step - 1])
        {
                #print "exit step - no improvement : $step \n";
                $step = $maxStep;
                $noImprovement = 1;
        }
}#end while step

close(out);

if($noImprovement)
{
        $stepCluster = $step - 2;
}
else
{
        $stepCluster = $step - 1;
}

#new cluster asignment
for($i = 0; $i < $#data; $i++)
{
#  $self->update_status("$moduleName: finding clusters - this may take a while (new cluster assignment step $i of ".$#data.")") ;
        $asignCluster[$i][0] = -1;
        $asignCluster[$i][1] = $pearsonThreshold;
        for($j = 0; $j < $clusters[$stepCluster]; $j++)
        {
                $meanX  = 0;
                $meanY  = 0;
                $meanXY = 0;
                $stdevX = 0;
                $stdevY = 0;
                $r      = 0;
                for($k = 1; $k <= $timepoints; $k++)
                {
                        $meanX  += $centroids[0][$i][$k];
                        $stdevX += $centroids[0][$i][$k] * $centroids[0][$i][$k];
                        
                        $meanY  += $centroids[$stepCluster][$j][$k];
                        $stdevY += $centroids[$stepCluster][$j][$k] * $centroids[$stepCluster][$j][$k];
                        
                        $meanXY += $centroids[0][$i][$k] * $centroids[$stepCluster][$j][$k];
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
                
                $r = ($meanXY - $meanX * $meanY) / ($stdevX * $stdevY);
                if ($r > $cluster[$i][1] && $i != $j)
                {
                        $asignCluster[$i][0] = $j;
                        $asignCluster[$i][1] = $r;
                }
        }
}

#print "Firepat: HC : ended ok \n";

return 1;
}
return 1;