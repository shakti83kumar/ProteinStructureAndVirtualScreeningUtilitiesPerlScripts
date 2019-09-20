#!usr/bin/perl -w
use strict;
use warnings;
# *************************************************
# * Copyright 2013 Shakti Kumar
# * 
# * CastpCavCrdExtraction.pl is free software: you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation, either version 3 of the License, or
# * (at your option) any later version.
# * 
# * LB-ADV_based_VS is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# * 
# * You should have received a copy of the GNU General Public License.
# * If not, see <http://www.gnu.org/licenses/>.
# *
# ***************************************************

my($CavityFileName, $RecptFolderName, $LigandFileName);
my(@PDBsNameArray, @LigFileFormatArray, $LigFlag, $RecptFlag);
my($YesOptionForLig, $LigPdbqtFolder, $YesOptionForPdb, $PdbToPdbqtFolder);

my $ArgumentNum = @ARGV;
#print("new argument: $ArgumentNum\n");
########### *docking Parameters* ###############
my $pdb2pdbqt_path = "/usr/share/pyshared/AutoDockTools/Utilities24";
#my @grid_dimension_array = (80, 80, 80);
my $num_modes = 5;
################################################
if($ArgumentNum == 6)
 {
   if(($ARGV[0] eq "-cfile")&&($ARGV[2] eq "-rdir")&&($ARGV[4] eq "-lfile"))
    {
       $CavityFileName  = $ARGV[1];
       $RecptFolderName = $ARGV[3];
       $LigandFileName  = $ARGV[5];    
       opendir(RECPDIR, "./$RecptFolderName")||die("\n\nUnable to open $RecptFolderName\nSome error in commond line input for receptor folder\nPlease see thorowly.............\nExiting the program.\n\n");
       @PDBsNameArray = grep(/\.pdb$/,readdir(RECPDIR));
       open(LIGFILE, "./$LigandFileName")||die("\n\nUnable to open $LigandFileName\nSome error in commond line input for ligand file.........\nPlease see thorowly.............\nExiting the program.\n\n");
       if(-z "./$LigandFileName")
         {
           die("Opss... $LigandFileName is empty file\n");
         }
       @LigFileFormatArray = split(/\./, $LigandFileName);
       $LigFlag   = &LigCheck(@LigFileFormatArray);
       $RecptFlag = &RecptCheck(@PDBsNameArray);
       print("ligand flag: $LigFlag and receptor flag: $RecptFlag\n");
       if(($LigFlag == 1)&&($RecptFlag == 1))
         {
           print("ligands are in correct file format and receptor folder contains pdb files\n");
           print("\n\n ***conversion of ligand file format to pdbqt format*** \n\n");
           print("Have you already converted all ligands in pdbqt format [if you have Enter 'YES'; else 'NO']:");
           $YesOptionForLig = <STDIN>;
           chop($YesOptionForLig);
           if($YesOptionForLig eq "YES")
             {
                print("Enter the pdbqt converted ligands containing folder name:");
                $LigPdbqtFolder = <STDIN>;
                chop($LigPdbqtFolder);
             }
           if($YesOptionForLig eq "NO")
             {
                $LigPdbqtFolder = &Lig2PdbqtPreparation($LigandFileName, $pdb2pdbqt_path);
             }
           print("\n\n ***conversion of receptor file format to pdbqt format*** \n\n"); 
           print("Have you already converted all proteins in pdbqt format [if you have Enter 'YES'; else 'NO']:");
           $YesOptionForPdb = <STDIN>;
           chop($YesOptionForPdb);
           if($YesOptionForPdb eq "YES")
             {
                print("Enter the pdbqt converted proteins containing folder name:");
                $PdbToPdbqtFolder = <STDIN>;
                chop($PdbToPdbqtFolder);
             }
           if($YesOptionForPdb eq "NO")
             {
                &Pdb2PdbqtPreparation($RecptFolderName, $pdb2pdbqt_path, @PDBsNameArray);
             }
           &ReadingCavityFile($CavityFileName, $RecptFolderName, $LigPdbqtFolder); 
          }
        if(($LigFlag == 0)&&($RecptFlag == 1))
          {
            die("\nyou have entered an unsupported file format of lignad\n",
                "file format of ligand should be either in mol2 or sdf\n\n");
          }
        if(($LigFlag == 1)&&($RecptFlag == 0))
          {
            die("Either file format of proteins is not in pdb format in receptor folders\nor No any pdb files are present in the folders\n");
          }
        if(($LigFlag == 0)&&($RecptFlag == 0))
          {
            die("niether file format of ligands are correct nor file format of receptor folder has pdb files\n");
          }
        
  
    } 
   else
    {
      &Help;
    }
close(LIGFILE);
closedir(RECPDIR);
 }
else
 {
   &Help;
 }

######################## *ligand checking* ############################
sub LigCheck
{
  my @LigandFileFormatArray = @_;
  my $LigandFlag = 0;
  if($LigandFileFormatArray[1] =~ /sdf|mol2/)
     {
        $LigandFlag = 1;
        return($LigandFlag);
     }
     else
     {
       return($LigandFlag);
     }
  unlink(@LigandFileFormatArray);
}

######################## *receptor checking* ##########################
sub RecptCheck
{
  my @ReceptNameArray = @_;
  my $RecptFlag = 0;
  if(@ReceptNameArray > 0)
     {
        $RecptFlag = 1;
        return($RecptFlag);
     }
  else 
     {
        return($RecptFlag);
     }
  unlink(@ReceptNameArray)
}

####################### *preparation of ligand files* ##################
sub Lig2PdbqtPreparation
{
    my ($ligfile, $prepare_lig_path) = @_;
    my ($model_line, $ligand);
    print("$ligfile\n");
    my @mol2file_array = split(/\./, $ligfile);
    my $smifile = $mol2file_array[0].".smi";
    my $pdbfile = $mol2file_array[0].".pdb";
    my $inputfile_format = $mol2file_array[1];
    print("$smifile, $pdbfile, $inputfile_format\n");
    print("converting the $inputfile_format to smile format\n");
    system("babel", "-i$inputfile_format", $ligfile, "-osmi", $smifile);
    print("converting the $inputfile_format to pdb format\n");
    system("babel", "-i$inputfile_format", $ligfile, "-opdb", $pdbfile);
    my $pdbfolder = $mol2file_array[0]."_splited_pdb";
    mkdir($pdbfolder, 0755);
    my $pdbqfolder = $mol2file_array[0]."_splited_pdbq";
    mkdir($pdbqfolder, 0755);
    my $pdbqtfolder = $mol2file_array[0]."_splited_pdbqt";
    mkdir($pdbqtfolder, 0755);
    open(LIGFILE, "$smifile");
    my $num = 1;
    while(my $smiline = <LIGFILE>)
        {
          $smiline =~ s/^\s+|\s+$//g;
          my @smiline_array = split(/\s+/, $smiline);
          my $ligname = $smiline_array[1];
          my $ligname_pdbfile = $smiline_array[1].".pdb";
          open(PDBOUTFILE, ">$pdbfolder/$ligname_pdbfile");
          open(PDBINFILE, "$pdbfile");
          my $flag = 0;
          while(my $mol2line = <PDBINFILE>)
               {
                 if($mol2line =~ /^MODEL/g)
                   {
                     $model_line = $mol2line;
                   }
                 if($mol2line =~ /$ligname/g)
                   {
                     print PDBOUTFILE ($model_line);
                     print("lig name: $ligname\n");
                     print("lig num: $num\n");
                     $flag = 1;
                   }
                 if($flag == 1)
                   { 
                     if($mol2line =~ /^ENDMDL/g)
                       {
                          last;
                       }
                     else
                       {
                           print PDBOUTFILE ($mol2line);
                       }
                   }
               }
          $num = $num + 1;
          close(PDBOUTFILE);
          close(PDBINFILE); 
        }
       opendir(PDBDIR, "./$pdbfolder");
       my @all_ligands = grep(/\.pdb$/, readdir(PDBDIR));
       foreach $ligand (@all_ligands)
             {
               my @lig_name_array = split(/\./, $ligand);
               my $lig_base_name = $lig_name_array[0];
               my $lig_base_name_q = $lig_name_array[0].".pdbq";
               my $lig_base_name_qt = $lig_name_array[0].".pdbqt";
               system("python", "$prepare_lig_path/prepare_ligand4.py", "-l", "./$pdbfolder/$ligand", "-o", "./$pdbqfolder/$lig_base_name_q");
               system("python", "$prepare_lig_path/pdbq_to_pdbqt.py", "-s", "./$pdbqfolder/$lig_base_name", "-o", "./$pdbqtfolder/$lig_base_name_qt");
             }
return($pdbqtfolder);
close(LIGFILE);
closedir(PDBDIR);
}

##################### *Preparation of proteins files* ##############################
sub Pdb2PdbqtPreparation
{
  my ($receptor_folder, $prepare_receptor_path, @PDBsNameArray) = @_;
  foreach my $pdbfile (@PDBsNameArray)
        {
           my @pdbfile_array = split(/\./, $pdbfile);
           my $pdbfile_basename = $pdbfile_array[0];
           my $prepared_receptor_name = $pdbfile_basename.".pdbqt";
           system("python", "$prepare_receptor_path/prepare_receptor4.py", "-r", "$receptor_folder/$pdbfile", "-o", "$receptor_folder/$prepared_receptor_name", "-A", "checkhydrogens");
        }
}

################### *Reading the Cavtity file* #####################################
sub ReadingCavityFile
{
      my($CavityFileName, $PdbToPdbqtFolder, $LigPdbqtFolder) = @_;
      open(CAVITYFILE, "./$CavityFileName");
      #open(OUTFILE, ">$InvlvedAtmInCavity");
      my $AtmFg   = 0;
      my $CumZcrd = 0;
      my $CumYcrd = 0;
      my $CumXcrd = 0;
      my $ENDflag = 0;
      my($line, $AtomID, $AtomType);
      my($AvgCrdX, $AvgCrdY, $AvgCrdZ);
      my(@pdbNameArray, @CavtyLineArray);
      my($ResName, $ChainID, $pdbline, $AtomName, $RsSeqNm, @RsSeqNmAry);
      my($ResdName, $PdbCnID, $ResSqNm, $Xcrd, $Ycrd, $Zcrd, $gamma, $eachelement);
      my($pdbName, $ProteinPdbqtFile, $GridDim, @XYZ_Crd_Array, @ResSqNmArray);
        
      while($line = <CAVITYFILE>)
            {
               $line =~ s/^\s+|\s+$//g;
               if(($line =~ /^\#/)&&($line =~ /pdbName/g))
                 {
                    @pdbNameArray = split(/=/,  $line);
                    $pdbName = $pdbNameArray[1];
                    print("\$pdbName:$pdbName\n");
                 }
               if(($line !~ /^\#/)&&($line !~ /^END/))
                 {    
                    @CavtyLineArray = split(/\s+/, $line);
                    $CavtyLineArray[1] =~ s/^\s+|\s+$//g;
                    $CavtyLineArray[2] =~ s/^\s+|\s+$//g;
                    $CavtyLineArray[3] =~ s/^\s+|\s+$//g;
                    $CavtyLineArray[4] =~ s/^\s+|\s+$//g;

                    $AtomID   = $CavtyLineArray[1];
                    $AtomType = $CavtyLineArray[2];
                    $ResName  = $CavtyLineArray[3];
                    $ChainID  = $CavtyLineArray[4];
                    #print("\$AtomID: $AtomID, \$AtomType: $AtomType, \$ResName:$ResName, \$ChainID:$ChainID\n");
                    #opendir(PDBDIR, "./$PdbToPdbqtFolder");
                    my $pdbNameFile = "./".$PdbToPdbqtFolder."/".$pdbName;
                    open(PDBNAME, "$pdbNameFile");
                    while($pdbline = <PDBNAME>)
                         {
                            $pdbline =~ s/^\s+|\s+$//g;
                            if($pdbline =~ /^ATOM/)
                              {
                                 $AtomName = substr($pdbline, 12, 4);
                                 $ResdName = substr($pdbline, 17, 3);
                                 $PdbCnID  = substr($pdbline, 21, 1);
                                 $ResSqNm  = substr($pdbline, 23, 4);
                                 $Xcrd     = substr($pdbline, 30, 8);
                                 $Ycrd     = substr($pdbline, 38, 8);
                                 $Zcrd     = substr($pdbline, 46, 8);
                 
                                 $ResSqNm  =~ s/^\s+|\s+$//g;
                                 $PdbCnID  =~ s/^\s+|\s+$//g;
                                 $ResdName =~ s/^\s+|\s+$//g;
                                 $AtomName =~ s/^\s+|\s+$//g;
                                 $Xcrd     =~ s/^\s+|\s+$//g;
                                 $Ycrd     =~ s/^\s+|\s+$//g;
                                 $Zcrd     =~ s/^\s+|\s+$//g;
                                 @ResSqNmArray = split("", $ResSqNm);
                                 #print("\@ResSqNmArray: @ResSqNmArray\n");
                                 $gamma = 0;
                                 foreach $eachelement (@ResSqNmArray)
                                   {
                                     if($eachelement =~ /\d/)
                                      {
                                          $RsSeqNmAry[$gamma] = $eachelement;
                                          $gamma++;
                                      }
                                   }
                                 $RsSeqNm = join("", @RsSeqNmAry);
                                 #print("\$RsSeqNm:$RsSeqNm\n");
                                 #print("$AtomName, $ResdName, $PdbCnID, $ResSqNm, $Xcrd, $Ycrd, $Zcrd\n");
                                 if(($AtomID == $RsSeqNm)&&($AtomType eq $AtomName)&&($ResName eq $ResdName)&&($ChainID eq $PdbCnID))
                                   {
                                      $CumXcrd = $CumXcrd + $Xcrd;
                                      $CumYcrd = $CumYcrd + $Ycrd;
                                      $CumZcrd = $CumZcrd + $Zcrd;
                                      $XYZ_Crd_Array[$AtmFg][0] = $Xcrd;
                                      $XYZ_Crd_Array[$AtmFg][1] = $Ycrd;
                                      $XYZ_Crd_Array[$AtmFg][2] = $Zcrd;
                                      $AtmFg++;
                                      #print("$AtmFg, $AtomID, $AtomType, $ResName, $ChainID\n");
                                      last;
                                   }
                              }
                         }
                 }  
               if($line =~ /^END/)
                 {
                    $ENDflag = $ENDflag + 1; 
                    $AvgCrdX = $CumXcrd/$AtmFg;
                    $AvgCrdY = $CumYcrd/$AtmFg;
                    $AvgCrdZ = $CumZcrd/$AtmFg;
                    $AvgCrdX = sprintf("%6.3f", $AvgCrdX);
                    $AvgCrdY = sprintf("%6.3f", $AvgCrdY);
                    $AvgCrdZ = sprintf("%6.3f", $AvgCrdZ);
                    $GridDim = &GridSizeCalculation($AtmFg, $AvgCrdX, $AvgCrdY, $AvgCrdZ, @XYZ_Crd_Array);
                    print("\$AvgCrdX: $AvgCrdX, \$AvgCrdY: $AvgCrdY, \$AvgCrdZ: $AvgCrdZ, \$GridDim: $GridDim\n");
                    my @pdbNameFileArray = split(/\./, $pdbName);
                    my $ProteinPdbqtFile = $pdbNameFileArray[0].".pdbqt";
                    my $ConfigFileName =  $pdbNameFileArray[0]."_$ENDflag".".conf";
                    open(CONFILE, ">$PdbToPdbqtFolder/$ConfigFileName");
                    print CONFILE ("receptor = $ProteinPdbqtFile\n");
                    print CONFILE ("center_x = $AvgCrdX\n");
                    print CONFILE ("center_y = $AvgCrdY\n");
                    print CONFILE ("center_z = $AvgCrdZ\n");
                    print CONFILE ("size_x = $GridDim\n");
                    print CONFILE ("size_y = $GridDim\n");
                    print CONFILE ("size_z = $GridDim\n");
                    print CONFILE ("num_modes = $num_modes\n");
                    &AutoDockVinaHTS($PdbToPdbqtFolder, $ProteinPdbqtFile, $LigPdbqtFolder, $ConfigFileName);
                    #@&MakeConfigFile($AvgCrdX, $AvgCrdY, $AvgCrdZ, $num_modes, @grid_dimension_array);
                    print("$AtmFg, $CumZcrd, $CumYcrd, $CumXcrd\n");
                    $AtmFg = 0;
                    $CumZcrd = $CumYcrd = $CumXcrd = 0;
                    @XYZ_Crd_Array = ();
                 }
            }
close(CAVITYFILE);
close(PDBNAME);
close(CONFILE);
}
###################### *help subroutine* ###########################################
sub Help
        {
            print("\n\nFor running this program please type following commond\n");
            print("perl CAVS.pl -cfile <cavity file name> -rdir <receptor folder name> -lfile <ligand file name>\n\n");
        }


###################### *high throughput drug sreening by autoDock vina tool* ###########
sub AutoDockVinaHTS
    {
       my ($PdbToPdbqtFolder, $ProteinPdbqtFile, $LigPdbqtFolder, $ConfigFileName) = @_;
       my $num = 1;
       my(@ligand_name_array, @ProteinPdbqtFileArray, @ConfigFileNameArray);
       my($DockingResultFolder);
       $DockingResultFolder = "DockingResult";
       mkdir($DockingResultFolder, 0755);
       @ProteinPdbqtFileArray = split(/\./, $ProteinPdbqtFile);
       @ConfigFileNameArray = split(/\./, $ConfigFileName);
       my $vina_docked_result_folder = $DockingResultFolder."/".$ProteinPdbqtFileArray[0]."_Vinadock_result";
       mkdir($vina_docked_result_folder, 0755);
       my $recpt_dock = $vina_docked_result_folder."/".$ConfigFileNameArray[0];
       mkdir($recpt_dock, 0755);
       my $Receptor = "./".$PdbToPdbqtFolder."/".$ProteinPdbqtFile;
       print("\$Receptor: $Receptor\n");
       $ConfigFileName = "./".$PdbToPdbqtFolder."/".$ConfigFileName;
       opendir(LIGDIR, "./$LigPdbqtFolder");
       @ligand_name_array = grep(/\.pdbqt$/, readdir(LIGDIR));
       foreach my $ligand_pdbqt_file (@ligand_name_array)
             {
                my @ligand_pdbqt_file_array = split(/\./, $ligand_pdbqt_file);
                my $ligand_file = "./".$LigPdbqtFolder."/".$ligand_pdbqt_file;
                my $ligand_base_name = $ligand_pdbqt_file_array[0];
                my $vina_output_result_pdbbqt_Name = $ConfigFileNameArray[0]."_".$ligand_base_name.".pdbqt";
                my $vina_output_result_path = $recpt_dock."/".$vina_output_result_pdbbqt_Name;
                my $vina_output_result_log_Name = $ConfigFileNameArray[0]."_".$ligand_base_name."_"."log";
                my $vina_output_log_path = $recpt_dock."/".$vina_output_result_log_Name;
                system("vina", "--receptor", $Receptor, "--ligand", $ligand_file, "--config", $ConfigFileName, "--out", $vina_output_result_path, "--log", $vina_output_log_path);
                print("ligand name is:$ligand_base_name\n");
                print("ligand nummber is:$num\n");
                print("#################################################################\n");
                print("#                                                               #\n");
                print("#  Thanks to SHAKTI KUMAR.....                                  #\n");
                print("#  for writting perl script for high throughput screening       #\n");
                print("#  by AutoDock Vina, CASTp and OpenBabel                        #\n");
                print("#                                                               #\n");
                print("#################################################################\n");
                $num = $num + 1;
             }
    }
closedir(LIGDIR);  
############################## *Grid Dimension Calcualtion* #####################################################
sub GridSizeCalculation
  {
    my($AtmLnFg, $AvgCrX, $AvgCrY, $AvgCrZ, @XYZ_Cr_Array) = @_;
    my($alpha, @SqdArray, $Sqd, @SortSqdArray);
    for($alpha = 0; $alpha < $AtmLnFg; $alpha++)
      {
        $Sqd = ($XYZ_Cr_Array[$alpha][0]-$AvgCrX)**2+($XYZ_Cr_Array[$alpha][1]-$AvgCrY)**2+($XYZ_Cr_Array[$alpha][2]-$AvgCrZ)**2;
        $Sqd = ($Sqd)**0.5;
        $SqdArray[$alpha] = $Sqd+0.5;
      }
    @SortSqdArray = sort{$b<=>$a}@SqdArray;
    #print("@SortSqdArray\n");
    return($SortSqdArray[0]);
    undef(@SortSqdArray);
    undef(@SqdArray);
  }
 

