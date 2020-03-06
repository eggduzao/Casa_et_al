
# Import
import os
import sys

#################################################
# BIN / STD
#################################################

# Region List
regionFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/0_Regions/Selected_Regions2.txt"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/1_Bin_Matrices/selected_regions_2/"
tl = "/scratch/egadegu/RCS/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/2_Subtraction_Matrices/selected_regions_2/"

# Open file
inputFileName = fl + "2_rcs.txt"
inFile = open(inputFileName, "w")

# Region Loop
regionFile = open(regionFileName, "rU")
for line in regionFile:

  ll = line.strip().split("\t")

  # Hic Name List
  cond1List = ["STAG1_P", "STAG2_P"]
  cond2List = ["STAG1_M", "STAG2_M"]
  outList = ["STAG1", "STAG2"]

  # Hic Name Loop
  for i in range(0,len(cond1List)):

    # Name
    hicname1 = "_".join([cond1List[i], ll[0], ll[1], ll[2], "bin"])
    hicname2 = "_".join([cond2List[i], ll[0], ll[1], ll[2], "bin"])
    outName1 = "_".join([outList[i], ll[0], ll[1], ll[2], "bin"])

    hicname3 = "_".join([cond1List[i], ll[0], ll[1], ll[2], "raw"])
    hicname4 = "_".join([cond2List[i], ll[0], ll[1], ll[2], "raw"])
    outName2 = "_".join([outList[i], ll[0], ll[1], ll[2], "raw"])

    # Condition List
    condList = ["25K_norm"]

    # Condition Loop
    for cond in condList:
 
      # Parameters
      if("25K" in cond): res = 25000
      elif("100K" in cond): res = 100000

      # Percentile List
      percList = ["80", "90"]

      # Percentile Loop
      for perc in percList:

        # Input
        resolution = str(res)
        inputMatrix1FileName = il + perc + "/" + cond + "/" + hicname1 + ".txt"
        inputMatrix2FileName = il + perc + "/" + cond + "/" + hicname2 + ".txt"
        tempLocation = tl + perc + "/" + cond + "/" + hicname1 + "/"
        outputMatrixPrefix = ol + perc + "/" + cond + "/" + outName1

        # Write
        inFile.write(" ".join([resolution, inputMatrix1FileName, inputMatrix2FileName, tempLocation, outputMatrixPrefix])+"\n")

        # Input
        resolution = str(res)
        inputMatrix1FileName = il + perc + "/" + cond + "/" + hicname3 + ".txt"
        inputMatrix2FileName = il + perc + "/" + cond + "/" + hicname4 + ".txt"
        tempLocation = tl + perc + "/" + cond + "/" + hicname3 + "/"
        outputMatrixPrefix = ol + perc + "/" + cond + "/" + outName2

        # Write
        inFile.write(" ".join([resolution, inputMatrix1FileName, inputMatrix2FileName, tempLocation, outputMatrixPrefix])+"\n")

inFile.close()


