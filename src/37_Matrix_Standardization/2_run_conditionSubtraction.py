
# Import
import os
import sys

"""
# Condition List
counter = 1
fl = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/6_Matrix_Standardization/input/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/6_Matrix_Standardization/1_matrices/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/6_Matrix_Standardization/2_subtractions/"
cond1List = ["STAG1_minus_auxin_binarized", "STAG1_minus_auxin_binarized", "STAG2_minus_auxin_binarized", "STAG1_minus_auxin_standardized", "STAG1_minus_auxin_standardized", "STAG2_minus_auxin_standardized"]
cond2List = ["STAG2_minus_auxin_binarized", "STAG1_plus_auxin_binarized", "STAG2_plus_auxin_binarized", "STAG2_minus_auxin_standardized", "STAG1_plus_auxin_standardized", "STAG2_plus_auxin_standardized"]

# Condition Loop
for i in range(0,len(cond1List)):

  # Name
  cond1 = cond1List[i]
  cond2 = cond2List[i]

  # Input
  resolution = "25000"
  inputMatrix1FileName = il + cond1 + ".txt"
  inputMatrix2FileName = il + cond2 + ".txt"
  tempLocation = "/scratch/eduardo/SUB/"
  outputMatrixPrefix = ol + cond1 + "--MINUS--" + cond2

  # Write
  inFile = open(fl + str(counter) + "_sub.txt", "w")
  inFile.write("\n".join([resolution, inputMatrix1FileName, inputMatrix2FileName, tempLocation, outputMatrixPrefix]))
  inFile.close()
  counter += 1
"""

# Condition List
il = "/home/egg/Projects/Wendt_Stag/Results/6_Matrix_Standardization/1_matrices/250K/"
ol = "/home/egg/Projects/Wendt_Stag/Results/6_Matrix_Standardization/2_subtractions/250K/"
cond1List = ["STAG1_plus_auxin_binarized", "STAG2_plus_auxin_binarized", "STAG1_plus_auxin_standardized", "STAG2_plus_auxin_standardized"]
cond2List = ["STAG1_minus_auxin_binarized", "STAG2_minus_auxin_binarized", "STAG1_minus_auxin_standardized", "STAG2_minus_auxin_standardized"]

# Condition Loop
for i in range(0,len(cond1List)):

  # Name
  cond1 = cond1List[i]
  cond2 = cond2List[i]

  # Input
  resolution = "250000"
  inputMatrix1FileName = il + cond1 + ".txt"
  inputMatrix2FileName = il + cond2 + ".txt"
  tempLocation = "./TEMP/"
  outputMatrixPrefix = ol + cond1 + "__MINUS__" + cond2

  # Write
  command = "python 2_conditionSubtraction.py "+" ".join([resolution, inputMatrix1FileName, inputMatrix2FileName, tempLocation, outputMatrixPrefix])
  os.system(command)


