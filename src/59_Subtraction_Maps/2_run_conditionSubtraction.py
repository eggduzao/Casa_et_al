
# Import
import os
import sys

#################################################
# BIN / STD
#################################################

# Hic Name List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/1_Bin_Std_Matrices/"
tl = "/scratch/egadegu/RSM/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/2_Subtraction_Matrices/"
cond1List = ["3B9_5plus_bin", "3B9_5plus_std", "69_127plus_bin", "69_127plus_std"]
cond2List = ["3B9_5-_bin", "3B9_5-_std", "69_127-_bin", "69_127-_std"]

# Open file
inputFileName = fl + "2_rsm.txt"
inFile = open(inputFileName, "w")

# Hic Name Loop
for i in range(0,len(cond1List)):

  # Name
  hicname1 = cond1List[i]
  hicname2 = cond2List[i]

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
      inputMatrix1FileName = il + perc + "/" + cond + "/" + hicname1 + ".txt"
      inputMatrix2FileName = il + perc + "/" + cond + "/" + hicname2 + ".txt"
      tempLocation = tl + perc + "/" + cond + "/" + hicname1 + "/"
      outputMatrixPrefix = ol + perc + "/" + cond + "/" + hicname1 + "_M_" + hicname2

      # Write
      inFile.write(" ".join([resolution, inputMatrix1FileName, inputMatrix2FileName, tempLocation, outputMatrixPrefix])+"\n")

#################################################
# RAW
#################################################

il = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/"
tl = "/scratch/egadegu/RSMR/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/2_Subtraction_Matrices/"
cond1List = ["3B9_5plus", "69_127plus"]
cond2List = ["3B9_5-", "69_127-"]

# Hic Name Loop
for i in range(0,len(cond1List)):

  # Name
  hicname1 = cond1List[i]
  hicname2 = cond2List[i]

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

    # Input
    resolution = str(res)
    inputMatrix1FileName = il + cond + "/" + hicname1 + ".txt"
    inputMatrix2FileName = il + cond + "/" + hicname2 + ".txt"
    tempLocation = tl + cond + "/" + hicname1 + "/"
    outputMatrixPrefix = ol + "raw/" + cond + "/" + hicname1 + "_M_" + hicname2

    # Write
    inFile.write(" ".join([resolution, inputMatrix1FileName, inputMatrix2FileName, tempLocation, outputMatrixPrefix])+"\n")

inFile.close()


