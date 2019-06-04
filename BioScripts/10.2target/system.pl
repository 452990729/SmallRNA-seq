#!/usr/bin/env perl

use strict;

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub AppendToFile
{
    system("echo \"$_[1]\" >> $_[0]");
}

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub CreateDir
{
  my $dir = $_[0];
  if (!(-d $dir)) { system("mkdir $dir"); }
}

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub CopyFile
{
  if (FileExists($_[0]))
  {
      system("cp $_[0] $_[1]");
  }
}

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub DeleteFile
{
  if (FileExists($_[0]))
  {
      system("rm $_[0]");
  }
}

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub FileExists
{
  if (-e $_[0]) { return 1; }
  else { return 0; }
}

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub FileIsADirectory
{
  if (-d $_[0]) { return 1; }
  else { return 0; }
}

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub GetNumLinesInFile
{
  my $r = int(rand 1000000000);

  system("wc $_[0] > tmp.$r");

  open(TMP, "<tmp.$r");
  my $line = <TMP>;
  chop $line;

  $line =~ /([^\s]+)/;

  &DeleteFile("tmp.$r");

  return $1;
}

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub LinkFile
{
  if (!(-l $_[1]))
  {
      system("ln -s $_[0] $_[1]");
  }
}

#-------------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------------
sub MoveFile
{
    if (FileExists($_[0]))
    {
	system("mv $_[0] $_[1]");
    }
}

1
