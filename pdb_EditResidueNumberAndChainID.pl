#!/usr/bin/perl -w
use strict;
use warnings;
my ($Inputpdb, $Line, $ResidueNum, $AddingNum, $ResidReNum, $Outputpdb, $ChainID, $numspace);
if($#ARGV < 0)
 {
   print "provide valid PDB file as input\n";
   exit;
 }
else
 {
   $Inputpdb  = $ARGV[0];
   $AddingNum = $ARGV[1];
   $ChainID   = $ARGV[2];
   $Inputpdb  = &Clean($Inputpdb);
   $AddingNum = &Clean($AddingNum);
   $Outputpdb = $Inputpdb.'-edit'.'.pdb';
   print ("$Inputpdb\n");
   open(OUTFILE, ">$Outputpdb");
   open(INFILE, "./$Inputpdb");
   while ($Line = <INFILE>)
    {
      if($Line =~ /^ATOM|^HETATM|^TER/)
       {
         #print($Line);
         $ResidueNum = substr($Line, 22, 4);
         $ResidueNum = &Clean($ResidueNum);
         $ResidReNum = $AddingNum + $ResidueNum;
         print ("$ResidReNum\n");
         substr($Line, 21, 1) = $ChainID;
         $numspace = 4-length($ResidReNum);
         substr($Line, 22, 4) = " "x$numspace.$ResidReNum;
         print OUTFILE ($Line);
       }
      else
       {
         print OUTFILE ($Line);
       }
    }
 }

sub Clean
{
   my $Val = shift;
   $Val =~ s/^\s+|\s+$//;
   return($Val);
}
