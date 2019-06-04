#!/usr/bin/env perl

use strict;
use FindBin '$Bin';

require "$Bin/load_args.pl";
require "$Bin/format_number.pl";
require "$Bin/system.pl";

if ($ARGV[0] eq "--help")
{
  print STDOUT <DATA>;
  exit;
}

my $file_ref;
my $file = $ARGV[0];
if (length($file) < 1 or $file =~ /^-/) 
{
  $file_ref = \*STDIN;
}
else
{
  open(FILE, $file) or die("Could not open file '$file'.\n");
  $file_ref = \*FILE;
}

my %args = load_args(\@ARGV);

my $include_sequence = get_arg("s", 0, \%args);
my $print_all = get_arg("a", 0, \%args);

while(<$file_ref>)
{
  chop;

  my @row = split(/\t/);

  print "$row[0]\t";
  print length($row[1]);
  
  if ($include_sequence == 1)
  {
    print "\t$row[1]";
  }
  elsif ($print_all == 1)
  {
     for (my $i = 1; $i <= @row; $i++)
     {
	print "\t$row[$i]";
     }
  }

  print "\n";
}

__DATA__

/PUBLIC/source/RNA/smallRNA/pipeline/sRNA/tools_script/PITA/lib/stab2length.pl <file>

   Outputs the length of each sequence

   -s: Include the sequence in the output file

   -a: Prints additional fields if exist

