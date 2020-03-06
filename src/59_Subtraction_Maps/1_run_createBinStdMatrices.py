
# Import
import os
import sys

# Condition List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/"
tl = "/scratch/egadegu/CBS/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/1_Bin_Std_Matrices/"
hicNameList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

# Open file
inputFileName = fl + "1_cbs.txt"
inFile = open(inputFileName, "w")

# Condition Loop
for hicName in hicNameList:

  # Condition List
  condList = ["25K_none", "25K_norm", "50K_none", "50K_norm", "100K_none", "100K_norm", "250K_none", "250K_norm"]

  # Condition Loop
  for cond in condList:

    # Parameters
    if("10K" in cond): res = 10000
    elif("25K" in cond): res = 25000
    elif("50K" in cond): res = 50000
    elif("100K" in cond): res = 100000
    elif("250K" in cond): res = 250000
    elif("500K" in cond): res = 500000
    elif("1000K" in cond): res = 1000000

    # Percentile List
    percList = ["70", "80", "90", "95", "98"]

    # Percentile Loop
    for perc in percList:

      # Input
      resolution = str(res)
      percentileThreshold = perc
      percentileFileName = il + cond + "/" + hicName + "_percentiles.txt"
      inputMatrixFileName = il + cond + "/" + hicName + ".txt"
      tempLocation = tl + perc + "/" + cond + "/" + hicName + "/"
      outBinPrefix = ol + perc + "/" + cond + "/" + hicName + "_bin"
      outStdPrefix = ol + perc + "/" + cond + "/" + hicName + "_std"

      # Write
      inFile.write(" ".join([resolution, percentileThreshold, percentileFileName, inputMatrixFileName, tempLocation, outBinPrefix, outStdPrefix])+"\n")

inFile.close()


