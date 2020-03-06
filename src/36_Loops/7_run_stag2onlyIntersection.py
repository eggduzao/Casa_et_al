
# Import
import os
import sys

# Loop List
sl = "/home/egg/Projects/Wendt_Stag/Data/stag_bed_files/"
ll = "/home/egg/Projects/Wendt_Stag/Results/5_Loops/3_contact_files_filtered/"
ol = "/home/egg/Projects/Wendt_Stag/Results/5_Loops/7_Stag2_Loopograms/loop_intersect_stag2_only/"
loopList = ["STAG1_minus_auxin_150Kb_5Mb" , "STAG1_plus_auxin_150Kb_5Mb" , "STAG2_minus_auxin_150Kb_5Mb" , "STAG2_plus_auxin_150Kb_5Mb"]

# Loop Loop
for loop in loopList:

  # Paramters
  name = "_".join(loop.split("_")[:-2])

  # Input
  stagRegionFileName = sl + "STAG2_only.bed"
  loopFileName = ll + loop + ".txt"
  tempLocation = "./TEMP/"
  outputFileName = ol + name + ".bed"

  # Execution
  command = "python 7_stag2onlyIntersection.py "+" ".join([stagRegionFileName, loopFileName, tempLocation, outputFileName])
  os.system(command)


