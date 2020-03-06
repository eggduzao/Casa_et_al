
# Import
import os
import sys

# Cell list
il1 = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/input/10K_norm/"
il2 = "/home/egg/Projects/Wendt_Stag/Previous_Results/3_Loopograms/matrix_new/10K_norm/"
ol = "/home/egg/Projects/Wendt_Stag/Previous_Results/13_Distance_Boxplots/"
inList = ["R_STAG1_RM_STAG1_intersection", "R_STAG1_RM_STAG1_minusOnly", "R_STAG1_RM_STAG1_plusOnly", "69_127plus.txt:STAG1_regions", "69_127-.txt:STAG1_regions"]

# Cell loop
for inName in inList:

  # Parameters
  il = il1
  suff = ".bed"
  if("txt" in inName):
    il = il2
    suff = ".txt"

  # Input
  loopRegionsFileName = il + inName + suff
  outputFileName = ol + inName + ".txt"

  # Execution
  command = "python 1_distanceTable.py "+" ".join([loopRegionsFileName, outputFileName])
  os.system(command)


