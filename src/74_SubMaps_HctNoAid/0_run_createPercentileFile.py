
# Import
import os
import sys

# Open file
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/"
inputFileName = fl + "0_cpf.txt"
inFile = open(inputFileName, "w")

# Input
inputMatrixFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/39_Process_Rao_Matrix/Untreated_Synchronized_25000.txt"
outputPercFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/39_Process_Rao_Matrix/Untreated_Synchronized_25000_percentiles.txt"

# Write
inFile.write(" ".join([inputMatrixFileName, outputPercFileName])+"\n")
inFile.close()


