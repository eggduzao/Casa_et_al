
# Import
import os
import sys

# Stag dictionary
sDict = {"3B9_5-": "STAG2_minus_auxin", "3B9_5plus": "STAG2_plus_auxin", "69_127-": "STAG1_minus_auxin", "69_127plus": "STAG1_plus_auxin"}

# Region List
rl = "/home/egg/Projects/Wendt_Stag/Results/5_Loops/7_Stag2_Loopograms/loop_intersect_stag2_only/"
ml = "/home/egg/Projects/Wendt_Stag/Data/stag_matrix_files/25K_norm/"
ol = "/home/egg/Projects/Wendt_Stag/Results/5_Loops/7_Stag2_Loopograms/loopogram_tables/"
regionList = ["STAG1_minus_auxin", "STAG1_plus_auxin", "STAG2_minus_auxin", "STAG2_plus_auxin"]

# Region Loop
for region in regionList:

  # Matrix List
  matrixList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

  # Matrix Loop
  for matrix in matrixList:

    # Input
    loopBins = "10"
    resolution = "25000"
    regionFileName = rl + region + ".bed"
    matrixFileName = ml + matrix + ".txt"
    outputMatrixFileName = ol + "R_" + region + "__M_" + sDict[matrix] + ".txt"
    
    command = "python 8_loopogram.py "+" ".join([loopBins, resolution, regionFileName, matrixFileName, outputMatrixFileName])
    os.system(command)


