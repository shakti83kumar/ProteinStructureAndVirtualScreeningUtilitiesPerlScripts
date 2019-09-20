#!/usr/bin/perl -w
use strict;
use warnings;
if($#ARGV < 0)
 {
   print "input the appropriate file\n";
   exit;
 }
else
 {
   my $input_name = $ARGV[0];
   $input_name =~ s/^\s+|\s+$//;
   my @input_name_array = split(/\./, $input_name);
   my $output_name = $input_name_array[0].'_edit'.'.'.$input_name_array[scalar(@input_name_array)-1];
   print $output_name, "\n";
   open(INPUT, $input_name);
   open(OUTPUT,">$output_name");
   while(my $line = <INPUT>)
    {
      $line =~ s/^\s+|\s+$//;
      $line =~ s/\s+/\t/g;
      print OUTPUT $line, "\n";
    }
 }
