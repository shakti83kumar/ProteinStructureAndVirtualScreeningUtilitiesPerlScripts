#!/usr/bin/perl
use strict;
use warnings;
my($nucs, $aas, $aasposition, $nucsposition);
my(%gencode, @blosum45, @blosum50, @blosum62);
my($FastaAlignmentFile, $SeqSimMatOutfile);
#
#
# input and output files
$FastaAlignmentFile = "Longest_modeled_cyst_protease_str_chain_added_clanCA_for_Alignment_group-2_result_20141015.afasta";
$SeqSimMatOutfile   = "Longest_modeled_cyst_protease_str_chain_added_clanCA_for_Alignment_group-2_result_20141015-ssmresults.rms_rot";
#
#
# The standard genetic code table 
%gencode = 
    (TTT => "F", TTC => "F", TTA => "L", TTG => "L",
     TCT => "S", TCC => "S", TCA => "S", TCG => "S",
     TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
     TGT => "C", TGC => "C", TGA => "*", TGG => "W",
     
     CTT => "L", CTC => "L", CTA => "L", CTG => "L",
     CCT => "P", CCC => "P", CCA => "P", CCG => "P", 
     CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
     CGT => "R", CGC => "R", CGA => "R", CGG => "R", 
     
     ATT => "I", ATC => "I", ATA => "I", ATG => "M",
     ACT => "T", ACC => "T", ACA => "T", ACG => "T", 
     AAT => "N", AAC => "N", AAA => "K", AAG => "K",
     AGT => "S", AGC => "S", AGA => "R", AGG => "R",

     GTT => "V", GTC => "V", GTA => "V", GTG => "V", 
     GCT => "A", GCC => "A", GCA => "A", GCG => "A",
     GAT => "D", GAC => "D", GAA => "E", GAG => "E",
     GGT => "G", GGC => "G", GGA => "G", GGG => "G" );

# The BLOSUM45 amino acid substitution matrix
@blosum45 = 
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
  ( [  5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-2,-2, 0],  # A
    [ -2, 7, 0,-1,-3, 1, 0,-2, 0,-3,-2, 3,-1,-2,-2,-1,-1,-2,-1,-2],  # R
    [ -1, 0, 6, 2,-2, 0, 0, 0, 1,-2,-3, 0,-2,-2,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-1, 2, 7,-3, 0, 2,-1, 0,-4,-3, 0,-3,-4,-1, 0,-1,-4,-2,-3],  # D
    [ -1,-3,-2,-3,12,-3,-3,-3,-3,-3,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1],  # C
    [ -1, 1, 0, 0,-3, 6, 2,-2, 1,-2,-2, 1, 0,-4,-1, 0,-1,-2,-1,-3],  # Q
    [ -1, 0, 0, 2,-3, 2, 6,-2, 0,-3,-2, 1,-2,-3, 0, 0,-1,-3,-2,-3],  # E
    [  0,-2, 0,-1,-3,-2,-2, 7,-2,-4,-3,-2,-2,-3,-2, 0,-2,-2,-3,-3],  # G
    [ -2, 0, 1, 0,-3, 1, 0,-2,10,-3,-2,-1, 0,-2,-2,-1,-2,-3, 2,-3],  # H
    [ -1,-3,-2,-4,-3,-2,-3,-4,-3, 5, 2,-3, 2, 0,-2,-2,-1,-2, 0, 3],  # I
    [ -1,-2,-3,-3,-2,-2,-2,-3,-2, 2, 5,-3, 2, 1,-3,-3,-1,-2, 0, 1],  # L
    [ -1, 3, 0, 0,-3, 1, 1,-2,-1,-3,-3, 5,-1,-3,-1,-1,-1,-2,-1,-2],  # K
    [ -1,-1,-2,-3,-2, 0,-2,-2, 0, 2, 2,-1, 6, 0,-2,-2,-1,-2, 0, 1],  # M
    [ -2,-2,-2,-4,-2,-4,-3,-3,-2, 0, 1,-3, 0, 8,-3,-2,-1, 1, 3, 0],  # F
    [ -1,-2,-2,-1,-4,-1, 0,-2,-2,-2,-3,-1,-2,-3, 9,-1,-1,-3,-3,-3],  # P
    [  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-3,-1,-2,-2,-1, 4, 2,-4,-2,-1],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1, 2, 5,-3,-1, 0],  # T
    [ -2,-2,-4,-4,-5,-2,-3,-2,-3,-2,-2,-2,-2, 1,-3,-4,-3,15, 3,-3],  # W
    [ -2,-1,-2,-2,-3,-1,-2,-3, 2, 0, 0,-1, 0, 3,-3,-2,-1, 3, 8,-1],  # Y
    [  0,-2,-3,-3,-1,-3,-3,-3,-3, 3, 1,-2, 1, 0,-3,-1, 0,-3,-1, 5]   # V 
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
    );
# The BLOSUM50 amino acid substitution matrix
@blosum50 = 
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
  ( [  5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0],  # A
    [ -2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3],  # R
    [ -1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4],  # D
    [ -1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1],  # C
    [ -1, 1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3],  # Q
    [ -1, 0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3],  # E
    [  0,-3, 0,-1,-3,-2,-3, 8,-2,-4,-4,-2,-3,-4,-2, 0,-2,-3,-3,-4],  # G
    [ -2, 0, 1,-1,-3, 1, 0,-2,10,-4,-3, 0,-1,-1,-2,-1,-2,-3, 2,-4],  # H
    [ -1,-4,-3,-4,-2,-3,-4,-4,-4, 5, 2,-3, 2, 0,-3,-3,-1,-3,-1, 4],  # I
    [ -2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5,-3, 3, 1,-4,-3,-1,-2,-1, 1],  # L
    [ -1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6,-2,-4,-1, 0,-1,-3,-2,-3],  # K
    [ -1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7, 0,-3,-2,-1,-1, 0, 1],  # M
    [ -3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8,-4,-3,-2, 1, 4,-1],  # F
    [ -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3],  # P
    [  1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5, 2,-4,-2,-2],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5,-3,-2, 0],  # T
    [ -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15, 2,-3],  # W
    [ -2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8,-1],  # Y
    [  0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5]   # V
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
  );
# The BLOSUM62 amino acid substitution matrix
@blosum62 = 
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
  ( [  4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],  # A
    [ -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],  # R
    [ -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],  # N
    [ -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],  # D
    [  0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],  # C
    [ -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],  # Q
    [ -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],  # E
    [  0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],  # G
    [ -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],  # H
    [ -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],  # I
    [ -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],  # L
    [ -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],  # K
    [ -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],  # M
    [ -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],  # F
    [ -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],  # P
    [  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],  # S
    [  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],  # T
    [ -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],  # W
    [ -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],  # Y
    [  0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4]   # V
    #  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
    );
# The amino acids, and a hashmap from amino acid to its index 0-19
#
$nucs = "ACGT";
$aas  = "ARNDCQEGHILKMFPSTWYV";
#
my($Line, @AlignmentArray, @FastaLineArray, $i, $j);
my(%FastaNamePosArray, @AllLineArray, @AllLine);
$i = $j = 0;
open(ALIFILE, "$FastaAlignmentFile");
open(OUTFILE, ">$SeqSimMatOutfile");
while(@AllLine = <ALIFILE>)
{
  @AllLineArray = @AllLine;
}
push (@AllLineArray, ">");
foreach $Line (@AllLineArray)
{
  if($Line !~ /^\s+/)
   {
     $Line =~ s/^\s+|\s+$//g;
     $AlignmentArray[$j] = $Line;
     if($Line =~ /^>/)
      {
        $Line =~ s/^>//g;
        $FastaLineArray[$i] = $Line;
        $FastaNamePosArray{$Line} = $j; 
        $i = $i + 1;
      }
    $j = $j + 1;
   }
}
undef($i);
undef($j);
undef($Line);
###########################***PFTE calculation***#######################################
my($lowlimtSeq1, $uplimitSeq1, $LengthSeq1, $Seq1, $k);
my($lowlimtSeq2, $uplimitSeq2, $LengthSeq2, $Seq2, $l);
my($PFTE, $antiPFTE, $p);
for($p = 0; $p<@FastaLineArray-1; $p++)
{
  print OUTFILE ("#  ", $p+1,"\t$FastaLineArray[$p]\n");
}
print OUTFILE ("####################*** Sequence Similarity Matrix ***##################\n");
for($k = 0; $k<@FastaLineArray-1; $k++)
{
  print         ($k+1,"\t$FastaLineArray[$k]\n");
  print OUTFILE ($k+1, "  ");
  $lowlimtSeq1 = $FastaNamePosArray{$FastaLineArray[$k]};
  $uplimitSeq1 = $FastaNamePosArray{$FastaLineArray[$k+1]};
  ($LengthSeq1, $Seq1) = &ExtractSeq($lowlimtSeq1, $uplimitSeq1, @AlignmentArray);
  for($l = 0; $l<@FastaLineArray-1; $l++)
    {
      $lowlimtSeq2 = $FastaNamePosArray{$FastaLineArray[$l]};
      $uplimitSeq2 = $FastaNamePosArray{$FastaLineArray[$l+1]};
      ($LengthSeq2, $Seq2) = &ExtractSeq($lowlimtSeq2, $uplimitSeq2, @AlignmentArray);
      $PFTE = SimilarityValue($LengthSeq1, $LengthSeq2, $Seq1, $Seq2, $aas);
      $PFTE = sprintf("%3.6f", $PFTE);
      $antiPFTE = 1.000000 - $PFTE;
      $antiPFTE = sprintf("%3.6f", $antiPFTE);
      print OUTFILE ("$antiPFTE  ");
      
    }
    print OUTFILE ("\n");
}
close(ALIFILE);
close(OUTFILE);
####################***To extract sequences***#######################
sub ExtractSeq
{
  my($low, $upp, @Alignment) = @_;
  my($Length, @OnlySeq, $Seq, @SeqArray, $EachChr);
  $Length = 0;
  $low = $low + 1;
  $upp = $upp - 1;
  @OnlySeq = @Alignment[$low..$upp];
  $Seq = join("", @OnlySeq);
  @SeqArray = split("", $Seq);
  foreach $EachChr (@SeqArray)
    {
      if($EachChr =~ /[A-Z]/)
       {
         $Length = $Length + 1;
       }
    }
    
  return($Length, $Seq);
}
####################***Smilarity value***##############################
sub SimilarityValue
{
  my($Length1, $Length2, $Sq1, $Sq2, $aar) = @_;
  my($LengthSq1, $LengthSq2, @Sq1Array, @Sq2Array);
  my($m, $Sq1Char, $Sq2Char, $aas1position, $pfte, $ssm);
  my($aas2position, $aaij, $NumAlignment, $aaii, $aajj);
  my($aaij2aaii, $aaij2aajj, $aaij2aaij);
  $LengthSq1 = length($Sq1);
  $LengthSq2 = length($Sq2);
  if($LengthSq1 == $LengthSq2)
   {
     $NumAlignment = $pfte = $ssm = 0;
     @Sq1Array = split("", $Sq1);
     @Sq2Array = split("", $Sq2);
     for($m = 0; $m<$LengthSq2; $m++)
       {
         if(($Sq1Array[$m] !~ /\-/)&&($Sq2Array[$m] !~ /\-/))
          {
            $Sq1Char = $Sq1Array[$m];
            $Sq2Char = $Sq2Array[$m];
            ($aas1position, $aas2position) = &CharPosition($Sq1Char, $Sq2Char, $aar);
            ##########***using of BLOSUM62 matrix***########## 
            $aaij = $blosum62[$aas1position][$aas2position];
            $aaii = $blosum62[$aas1position][$aas1position];
            $aajj = $blosum62[$aas2position][$aas2position];
            $aaij2aaii = $aaij/$aaii;
            $aaij2aajj = $aaij/$aajj;
            $aaij2aaij = ($aaij2aaii+$aaij2aajj)/2;
            $ssm = $ssm + ($aaij2aaij);
            $NumAlignment = $NumAlignment + 1;
            
          }
       }
   }
$ssm = $ssm/$NumAlignment;
return($ssm);
}
#######################*** charracter positions***##############################
sub CharPosition
{
  my($char1, $char2, $aastr) = @_;
  my($char1position, $char2position);
  my($aastr1, $aastr2);
  $aastr1 = $aastr2 = $aastr;
  $aastr1 =~ /$char1/g;
  $char1position = pos($aastr1);
  $char1position = $char1position - 1;
  $aastr2 =~ /$char2/g; 
  $char2position = pos($aastr2);
  $char2position = $char2position - 1;
  return($char1position, $char2position);
}
