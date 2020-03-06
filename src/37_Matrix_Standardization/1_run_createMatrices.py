
# Import
import os
import sys

"""
# Condition List
counter = 1
fl = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/6_Matrix_Standardization/input/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_matrix_files/25K_norm/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/6_Matrix_Standardization/1_matrices/"
condNameList = ["69_127-", "69_127plus", "3B9_5-", "3B9_5plus"]
condLabelList = ["STAG1_minus_auxin", "STAG1_plus_auxin", "STAG2_minus_auxin", "STAG2_plus_auxin"]

# Condition Loop
for i in range(0,len(condNameList)):

  # Name
  name = condNameList[i]
  label = condLabelList[i]

  # Input
  resolution = "25000"
  inputMatrixFileName = il + name + ".txt"
  tempLocation = "/scratch/eduardo/CBS/"
  outBinPrefix = ol + label + "_binarized"
  outStdPrefix = ol + label + "_standardized"

  # Write
  inFile = open(fl + str(counter) + "_cbs.txt", "w")
  inFile.write("\n".join([resolution, inputMatrixFileName, tempLocation, outBinPrefix, outStdPrefix]))
  inFile.close()
  counter += 1
"""

# Condition List
il = "/home/egg/Projects/Wendt_Stag/Data/stag_matrix_files/250K_norm/"
ol = "/home/egg/Projects/Wendt_Stag/Results/6_Matrix_Standardization/1_matrices/250K/"
condNameList = ["69_127-", "69_127plus", "3B9_5-", "3B9_5plus"]
condLabelList = ["STAG1_minus_auxin", "STAG1_plus_auxin", "STAG2_minus_auxin", "STAG2_plus_auxin"]

# Condition Loop
for i in range(0,len(condNameList)):

  # Name
  name = condNameList[i]
  label = condLabelList[i]

  # Input
  resolution = "250000"
  percentileThreshold = "25"
  percentileFileName = il + name + "_percentiles.txt"
  inputMatrixFileName = il + name + ".txt"
  tempLocation = "./TEMP/"
  outBinPrefix = ol + label + "_binarized"
  outStdPrefix = ol + label + "_standardized"

  # Write
  command = "python 1_createMatrices.py "+" ".join([resolution, percentileThreshold, percentileFileName, inputMatrixFileName, tempLocation, outBinPrefix, outStdPrefix])
  os.system(command)


