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
###
# split into components [html output]

package SrnaTools::Module::FirepatSplit;
use base SrnaTools::Module ;

#use strict;
#use warnings;

sub run{
my $self = shift ;

# FiRePat = find regulatory patterns
my $moduleName = "FiRePat - split";

$in     = $self->param('infile');
$ts     = $self->param('timeseries');
$outdir = $self->param('outdir');

$self->update_status("$moduleName: creating html output") ;   

open ts, $ts or SrnaTools::Exception::ModuleExecution->throw(message=>"cannot open timeseries");
open in, $in or SrnaTools::Exception::ModuleExecution->throw(message=>"cannot open correlation file");

while(<ts>)
{
        chomp;
        @d = split/,/;
        $timepoints = $#d;
}
close(ts);
open out1, '>',$outdir.'corr_in1.csv' or SrnaTools::Exception::ModuleExecution->throw(message=>"cannot open out1");
open out2, '>',$outdir.'corr_in2.csv' or SrnaTools::Exception::ModuleExecution->throw(message=>"cannot open out2");

while(<in>)
{
        chomp;
        @d = split/,/;
        
        for($i = 0; $i <= $timepoints; $i++)
        {
                print out1 "$d[$i],";
        }
        print out1 "\n";
        
        for($i = 0; $i <= $timepoints; $i++)
        {
                print out2 "$d[$i + $timepoints + 1],";
        }
        print out2 "\n";
}
close(in);
close(out1);
close(out2);

#print "Firepat:Split: script ended ok \n";

return 1;
}

return 1;