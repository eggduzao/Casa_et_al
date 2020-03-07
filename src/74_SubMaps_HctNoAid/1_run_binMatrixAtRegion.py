
# Import
import os
import sys

# Region List
regionFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/0_Regions/Selected_Regions.txt"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/"
il1 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/25K_norm/"
il2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/39_Process_Rao_Matrix/"
tl = "/usr/users/egadegu/scratch/CBF/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/1_Bin_Matrices/"

# Open file
inputFileName = fl + "1_cbf.txt"
inFile = open(inputFileName, "w")

# Region Loop
regionFile = open(regionFileName, "rU")
for line in regionFile:

  ll = line.strip().split("\t")

  # Hic List
  hicFileNameList = [il2 + "Untreated_Synchronized_25000", il1 + "69_127plus", il1 + "3B9_5plus"]
  newHicFileNameList = ["Untreated", "STAG1+", "STAG2+"]

  # Hic Loop
  for i in range(0,len(hicFileNameList)):

    # Parameters
    hicName = hicFileNameList[i]
    newHicName = newHicFileNameList[i]
    cond = "25K_norm"
    res = "25000"
    perc = "90"
    outName = "_".join([newHicName, ll[0], ll[1], ll[2]])

    # Input
    chromosome = ll[0]
    region1 = ll[1]
    region2 = ll[2]
    resolution = res
    percentileThreshold = perc
    percentileFileName = hicName + "_percentiles.txt"
    inputMatrixFileName = hicName + ".txt"
    tempLocation = tl + outName + "/"
    outBinPrefix = ol + outName + "_bin"
    outRawPrefix = ol + outName + "_raw"

    # Write
    inFile.write(" ".join([chromosome, region1, region2, resolution, percentileThreshold, percentileFileName, inputMatrixFileName, tempLocation, outBinPrefix, outRawPrefix])+"\n")

inFile.close()


