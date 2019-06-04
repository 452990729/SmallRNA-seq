#!/usr/bin/perl

###############################################################
#
# This script assist in adding downtime records into config/downtime_schedule.txt
# The site automatically checks that file for scheduled downtimes and
# does not accept jobs during those times.
# The script also deletes any records in the file that are in 
# the past.
# Format of time ranges in downtime_schedule.txt:
# DD-MM-YYYY;HH:MM DD-MM-YYYY;HH:MM
# where first time is the start and second the end of the range
# 
# Ranges can also overlap, so for example to extend a current downtime
# it is possible to simply add a new downtime from now to 
# some time in the future without deleting the currently active
# downtime record first.

use strict;
use warnings;
use Time::Piece;
use Time::Seconds;
use File::Basename;
use File::Spec;

#############################
# Declarations
#############################

# Assume that this script is in lib and that
# config is in top level dir of the site
my $schedule_file_dir = File::Spec->catdir(dirname($0), '..','config');
my $schedule_file = File::Spec->catfile($schedule_file_dir, "downtime_schedule.txt") ;

my $time_now = localtime;
my $time_now_plus_one = $time_now + ONE_DAY ;
my @downtime_records = () ;

#############################
# Check if the schedule file is
# present. If it is and it is not
# empty, read scheduled
# times in array if they are 
# not already in the past. This
# will later be printed back to 
# the file again
#############################

die "Could not find  directory '$schedule_file_dir'. Please create this directory first before running this script." if (! -d $schedule_file_dir) ;

if (-e $schedule_file && !-z $schedule_file) {
  # file exists and is not empty
  # Read contents and keep if end
  # of range still in the future
  open (SF, '<', $schedule_file) || die "Could not open $schedule_file\n" ;
  while (<SF>) {
    chomp ;
    my $end_time_str ;
    my $end_time ;
    if (/^.+? (.+)$/) {
      $end_time_str = $1 ;
    } else {
      die "Encountered incorrectly formatted line in schedule file: $_\n" ;
    }
    eval {$end_time = Time::Piece->strptime($end_time_str, "%d-%m-%Y;%H:%M")};
    die "Encountered incorrectly formatted line in schedule file: $end_time_str\n" if $@;
    push @downtime_records, $_ if $end_time >= $time_now ;
  }
  close SF ;
  
}



print <<HEREDOC ;

###############################################################
#                                                             
# This script adds an entry to the scheduled cluster downtime 
# file /config/downtime_schedule.txt. This file is read by    
# input_form.cgi, which will not accept jobs if the current   
# time overlaps with one of the scheduled downtimes.          
#                                                             
# It is ok to enter overlapping downtime ranges, e.g. to      
# extent a currently active downtime.                         
#                                                             
# The script asked for start and end time of the range.       
# You can hit return at each prompt to use the current        
# day, hour, minute etc. (from when the script was started)   
# The default for the end time is now + one day               
#                                                             
# Time on this machine now: $time_now                          
#                                                             
###############################################################

HEREDOC

#############################
# If we have any current or
# future down times, ask user
# if we should keep them
#############################

if (@downtime_records) {
  print "There are some current or future downtime records.\n" ;
  print "Do you want to confirm/delete them (otherwise all will be kept) [y/n]?" ;
  if (<> =~/y/i) {
    print "\nRecords are in the format 'DD-MM-YY;HH:MM DD-MM-YY;HH:MM' \nwhere the first time is the begin and the second the end time\n\n" ;
    my @temp_recs ;
    foreach (@downtime_records) {
      print "keep '$_' [y/n]?";
      if (<> =~/y/i) {
        print "  keeping record\n" ;
        push @temp_recs, $_ ;
      } else {
        print "  removing record\n" ;
      }
    }
    @downtime_records = @temp_recs ;
    
    print "Do you want to proceed to adding new records [y/n]?" ;
    unless  (<> =~/y/i) {
      write_to_file($schedule_file, \@downtime_records) ;
      exit 0 ;
    }
  }
}


#############################
# Ask user input.
#############################

my %time_range = (
  begin => {
    day   => $time_now->mday,
    month => $time_now->mon,
    year  => $time_now->yy,
    hour  => $time_now->hour,
    minute=> $time_now->min
  },
  end => {
    day   => $time_now_plus_one->mday,
    month => $time_now_plus_one->mon,
    year  => $time_now_plus_one->yy,
    hour  => $time_now_plus_one->hour,
    minute=> $time_now_plus_one->min
  },
) ;
my @sorted_keys= qw(day month year hour minute) ;
foreach my $tp (sort keys %time_range) {
  print "Enter $tp time:\n" ;
  my $is_valid = 0 ;
  while (!$is_valid) {
    my %in_time ;
    my @test_time ; # for validating format
    foreach my $param (@sorted_keys) {
      print "  $param [".$time_range{$tp}{$param}."] ?" ;
      chomp(my $in = <>) ; 
      $in_time{$param} = $in ne '' ? $in : $time_range{$tp}{$param};
      push @test_time, $in_time{$param} ;
    }
    my $time_str = join(' ', @test_time) ;
    eval {
      $time_range{$tp}{time_obj} = Time::Piece->strptime($time_str, "%d %m %Y %H %M");
    } ;
    # Check that the time could be parsed
    $is_valid = 1 unless $@ ;
    if ($is_valid) {
      print "  --> Time entered: $time_range{$tp}{time_obj}\n" ;
      # Check that end time isn't already in the past
      if ($tp eq 'end' &&  $time_range{$tp}{time_obj} < $time_now ) {
        $is_valid = 0 ;
        print ">>>> Time is valid but end time is in the past, please try again\n" ;
      }
      if ($tp eq 'end' && $time_range{end}{time_obj} <= $time_range{begin}{time_obj}) {
        $is_valid = 0 ;
        print ">>>> Time is valid but end time is before begin time, please try again\n" ;
      }
    } else {
      print ">>>> The format was invalid, please try again\n" ;
    }
  }
}

my @times ;
foreach my $tp (sort keys %time_range) {
  push @times, $time_range{$tp}{time_obj}->strftime("%d-%m-%y;%H:%M") ;
}
my $times_str = join(' ', @times) ;
print "\nScheduled new down time: $times_str\n" ;
push @downtime_records, $times_str ;
write_to_file($schedule_file, \@downtime_records) ;
exit 0 ;


###### subs

#############################
# write to file 
#############################
sub write_to_file{
  my $schedule_file = shift ;
  my $downtime_records_ref = shift ;
  open (SF, '>', $schedule_file) || die "Could not open $schedule_file\n" ;
  foreach (@$downtime_records_ref) {
    print SF $_,"\n" ;
  }
  close SF ;
  print "\nFile saved\n\n" ;
}

