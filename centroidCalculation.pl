#!/usr/bin/perl -w
use strict;
use warnings;
my $pdbname = $ARGV[0]; ## input PDB name
my($pdbline, $PdbCnID, $ResSqNm, $Xcrd, $Ycrd, $Zcrd);
my($AvgCrdX, $AvgCrdY, $AvgCrdZ); 
my $AtmFg   = 0;
my $CumZcrd = 0;
my $CumYcrd = 0;
my $CumXcrd = 0;
      open(PDBNAME, $pdbname)||die("Input file is not available\n");
      while($pdbline = <PDBNAME>)
        {
          $pdbline =~ s/^\s+|\s+$//g;
          if($pdbline =~ /^ATOM|^HETATM/)
           {
             $PdbCnID  = substr($pdbline, 21, 1);
             $ResSqNm  = substr($pdbline, 23, 4);
             $Xcrd     = substr($pdbline, 30, 8);
             $Ycrd     = substr($pdbline, 38, 8);
             $Zcrd     = substr($pdbline, 46, 8);
             $PdbCnID  =~ s/^\s+|\s+$//g;
             $ResSqNm  =~ s/^\s+|\s+$//g;
             $Xcrd     =~ s/^\s+|\s+$//g;
             $Ycrd     =~ s/^\s+|\s+$//g;
             $Zcrd     =~ s/^\s+|\s+$//g;
             $CumXcrd = $CumXcrd + $Xcrd;
             $CumYcrd = $CumYcrd + $Ycrd;
             $CumZcrd = $CumZcrd + $Zcrd;
             $AtmFg++;
             #print("$AtmFg, $AtomID, $AtomType, $ResName, $ChainID\n");
           }
        }
       $AvgCrdX = $CumXcrd/$AtmFg;
       $AvgCrdY = $CumYcrd/$AtmFg;
       $AvgCrdZ = $CumZcrd/$AtmFg;
       $AvgCrdX = sprintf("%6.3f", $AvgCrdX);
       $AvgCrdY = sprintf("%6.3f", $AvgCrdY);
       $AvgCrdZ = sprintf("%6.3f", $AvgCrdZ);
#       print OUTNAME ("\$AvgCrdX: $AvgCrdX, \$AvgCrdY: $AvgCrdY, \$AvgCrdZ: $AvgCrdZ\n");
       print ("\$AvgCrdX: $AvgCrdX, \$AvgCrdY: $AvgCrdY, \$AvgCrdZ: $AvgCrdZ\n");
       close(PDBNAME);
