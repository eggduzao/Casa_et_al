
# Import
import os
import sys

# Loop List
rl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_bed_files/"
ll = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/5_Loops/3_contact_files_filtered/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/5_Loops/5_ChoosingBestLoops/Venn/"
loopList = ["STAG1_minus_auxin_50Kb_5Mb", "STAG1_minus_auxin_50Kb_10Mb", "STAG1_minus_auxin_150Kb_5Mb" , "STAG1_minus_auxin_150Kb_10Mb",
            "STAG1_plus_auxin_50Kb_5Mb", "STAG1_plus_auxin_50Kb_10Mb", "STAG1_plus_auxin_150Kb_5Mb" , "STAG1_plus_auxin_150Kb_10Mb",
            "STAG2_minus_auxin_50Kb_5Mb", "STAG2_minus_auxin_50Kb_10Mb", "STAG2_minus_auxin_150Kb_5Mb" , "STAG2_minus_auxin_150Kb_10Mb",
            "STAG2_plus_auxin_50Kb_5Mb", "STAG2_plus_auxin_50Kb_10Mb", "STAG2_plus_auxin_150Kb_5Mb" , "STAG2_plus_auxin_150Kb_10Mb"]
regionList = ["STAG1_only", "STAG2_only", "shared"]

# Loop Loop
for loop in loopList:

  # Input
  resolution = "25000"
  loopFileName = ll + loop + ".txt"
  stagRegionList = ",".join([rl + e + ".bam" for e in regionList])
  tempLocation = ol + "TEMP/"
  outputFileName = ol + loop + ".txt"

  # Execution
  command = "python 6_vennDiagram.py "+" ".join([resolution, loopFileName, stagRegionList, tempLocation, outputFileName])
  os.system(command)


