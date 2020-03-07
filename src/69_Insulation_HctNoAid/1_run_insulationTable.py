
# Import
import os
import sys

# Region List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/40_Fig4G_Insulation_RaoNeg/input/"
rl = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bed_files/"
ml1 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/10K_norm/"
ml2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/39_Process_Rao_Matrix/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/40_Fig4G_Insulation_RaoNeg/1_Insulation_Table/"
regionList = ["STAG1_only", "STAG2_only", "shared"]

# Write
inFileName = fl + "1_cip.txt"
inFile = open(inFileName, "w")

# Region Loop
for regionName in regionList:

  # Matrix List
  if("STAG1" in regionName):
    matrixListList = [[ml2+"Untreated_10000"], [ml2+"Untreated_Synchronized_10000"], [ml1+"69_127plus", ml2+"Untreated_10000"], [ml1+"69_127plus", ml2+"Untreated_Synchronized_10000"]]
    labelList = ["Untreated", "Synchronized", "STAG1_DEG_minus_Untreated", "STAG1_DEG_minus_Synchronized"]
  elif("STAG2" in regionName):
    matrixListList = [[ml2+"Untreated_10000"], [ml2+"Untreated_Synchronized_10000"], [ml1+"3B9_5plus", ml2+"Untreated_10000"], [ml1+"3B9_5plus", ml2+"Untreated_Synchronized_10000"]]
    labelList = ["Untreated", "Synchronized", "STAG2_DEG_minus_Untreated", "STAG2_DEG_minus_Synchronized"]
  else:
    matrixListList = [[ml2+"Untreated_10000"], [ml2+"Untreated_Synchronized_10000"], [ml1+"69_127plus", ml2+"Untreated_10000"], [ml1+"69_127plus", ml2+"Untreated_Synchronized_10000"], [ml1+"3B9_5plus", ml2+"Untreated_10000"], [ml1+"3B9_5plus", ml2+"Untreated_Synchronized_10000"]]
    labelList = ["Untreated", "Synchronized", "STAG1_DEG_minus_Untreated", "STAG1_DEG_minus_Synchronized", "STAG2_DEG_minus_Untreated", "STAG2_DEG_minus_Synchronized"]

  # Stag Loop
  for i in range(0,len(matrixListList)):

    # Parameters
    matrixList = matrixListList[i]
    label = labelList[i]

    # Input
    halfBin = "20"
    resolution = "10000"
    regionFileName = rl + regionName + ".bed"
    matrix1FileName = matrixList[0] + ".txt"
    if(len(matrixList) > 1): matrix2FileName = matrixList[1] + ".txt"
    else: matrix2FileName = "NA"
    outputTableFileName = ol + regionName + "/" + label + ".txt"

    # Write
    inFile.write(" ".join([halfBin, resolution, regionFileName, matrix1FileName, matrix2FileName, outputTableFileName])+"\n")

inFile.close()


