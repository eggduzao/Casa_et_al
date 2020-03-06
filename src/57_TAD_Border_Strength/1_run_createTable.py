
# Import
import os
import sys

# TAD List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/28_TAD_Border_Strength/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/4_TADs/1_TADs_GMAP/"
ml = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/25K_norm/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/28_TAD_Border_Strength/"
tadList = ["STAG1-AUX", "STAG1+AUX", "STAG2-AUX", "STAG2+AUX"]
matList = ["69_127-", "69_127plus", "3B9_5-", "3B9_5plus"]

# Open File
inputFileName = fl + "1_ctb.txt"
inFile = open(inputFileName, "w")

# TAD Loop
for i in range(0,len(tadList)):

  # Parameter
  tadfile = tadList[i]
  matfile = matList[i]

  # Input
  resolution = "25000"
  matrixFileName = ml + matfile + ".txt"
  tadFileName = il + tadfile + ".txt"
  outputFileName = ol + tadfile + ".txt"

  # Execution
  inFile.write(" ".join([resolution, matrixFileName, tadFileName, outputFileName])+"\n")

inFile.close()


