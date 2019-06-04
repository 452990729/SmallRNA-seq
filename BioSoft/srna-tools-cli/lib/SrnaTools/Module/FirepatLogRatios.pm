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
# It calculates log2 ratios for FiRePat.
# TODO This class has been converted from a standalone script and might
# still need some cleanup 
###

package SrnaTools::Module::FirepatLogRatios;
use base SrnaTools::Module ;

#use strict;
#use warnings;

sub run{
my $self = shift ;

my $moduleName = "FiRePat - log2-ratios";
#compute log2 ratios relative to first timepoint
$in = $self->param('infile');
$out= $self->param('outfile');
$self->update_status("$moduleName: calculating log2 ratios") ;

open in, $in or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open input file file");
open out, '>', $out or SrnaTools::Exception::ModuleExecution->throw(message=>"$moduleName, can not open output file file");

$#data = -1;
while(<in>)
{       
        chomp;
        @d = split/,/;
        $tp = $#d;
        push @data,[@d];
}
close(in);

if($#data < 0)
{
	$self->job->_run_module(
		'firepat_create_failure_outfiles',
		{
			message => 'Firepat: CreateLogRatios:empty input file.'
		}
	);
	return 0;
}

$countOut = 0;
$eps = 0.001;
for($i = 0; $i <= $#data; $i++)
{
        $memo = $data[$i][1];
        if($memo == 0)
        {
                $memo = 1;
        }
        print out "$data[$i][0],";
        print out "0,";
		$countOut++;
        for($j = 2; $j <= $tp; $j++)
        {
                if($data[$i][$j-1] < $eps)
                {
                        $data[$i][$j-1] = $eps;
                }
                if($data[$i][$j] < 0.001)
                {
                   $value = log(($data[$i][$j] + 0.001)/ $data[$i][$j-1])/log(2);
                }
                else
                {
                   $value = log(($data[$i][$j])/ $data[$i][$j-1])/log(2);
                }
                print out "$value, ";
        }
        print out "\n";
}
close(out);

return 1;
}

return 1;