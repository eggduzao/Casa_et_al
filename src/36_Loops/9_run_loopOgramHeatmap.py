
# Import
import os
import sys

# Region List
tl = "/home/egg/Projects/Wendt_Stag/Results/5_Loops/7_Stag2_Loopograms/loopogram_tables/"
rl = "/home/egg/Projects/Wendt_Stag/Results/5_Loops/7_Stag2_Loopograms/loop_intersect_stag2_only/"
ol = "/home/egg/Projects/Wendt_Stag/Results/5_Loops/7_Stag2_Loopograms/loopogram_plots/"
regionList = ["STAG1_minus_auxin", "STAG1_plus_auxin", "STAG2_minus_auxin", "STAG2_plus_auxin"]

# Region Loop
for region in regionList:

  # Matrix Loop
  matrixList = ["STAG1_minus_auxin", "STAG1_plus_auxin", "STAG2_minus_auxin", "STAG2_plus_auxin"]

  # Region Loop
  for matrix in matrixList:

    # Input
    inputTableFileName = tl + "R_" + region + "__M_" + matrix + ".txt"
    regionsFileName = rl + region + ".bed"
    outputFileName = ol + "R_" + region + "__M_" + matrix + ".pdf"

    command = "R CMD BATCH '--args '"+inputTableFileName+"' '"+regionsFileName+"' '"+outputFileName+" 9_loopOgramHeatmap.R 9_loopOgramHeatmap.Rout"
    os.system(command)


