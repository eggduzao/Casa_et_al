
# Import
import os
import sys

#################################################
# BIN / STD
#################################################

# Region List
regionFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/0_Regions/Selected_Regions.txt"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/1_Bin_Matrices/"
tl = "/usr/users/egadegu/scratch/RCS/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/2_Subtraction_Matrices/"

# Open file
inputFileName = fl + "2_rcs.txt"
inFile = open(inputFileName, "w")

# Region Loop
regionFile = open(regionFileName, "rU")
for line in regionFile:

  ll = line.strip().split("\t")

  # Hic Name List
  cond1List = ["STAG1+", "STAG2+"]
  cond2List = ["Untreated", "Untreated"]
  outList = ["STAG1+_v_Untreated", "STAG2+_v_Untreated"]

  # Hic Name Loop
  for i in range(0,len(cond1List)):

    # Parameters
    hicname1 = "_".join([cond1List[i], ll[0], ll[1], ll[2], "raw"])
    hicname2 = "_".join([cond2List[i], ll[0], ll[1], ll[2], "raw"])
    outName1 = "_".join([outList[i], ll[0], ll[1], ll[2], "raw"])
    res = "25000"

    # Input
    resolution = res
    inputMatrix1FileName = il + hicname1 + ".txt"
    inputMatrix2FileName = il + hicname2 + ".txt"
    tempLocation = tl + outName1 + "/"
    outputMatrixPrefix = ol + outName1

    # Write
    inFile.write(" ".join([resolution, inputMatrix1FileName, inputMatrix2FileName, tempLocation, outputMatrixPrefix])+"\n")

inFile.close()


