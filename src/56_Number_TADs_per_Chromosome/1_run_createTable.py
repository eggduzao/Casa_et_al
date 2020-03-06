
# Import
import os
import sys

# Open file
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/27_Number_TADs_per_Chromosome/input/"
inputFileName = fl + "1_ctt.txt"
inFile = open(inputFileName, "w")

# Input
prefix = "T_3_80_25_100_5_10_0.95_0.5"
categoryList = "69_127-,69_127plus,3B9_5-,3B9_5plus"
inputLocation = "/usr/users/egadegu/Projects/Wendt_Stag/Results/4_TADs/1_TADs_GMAP/1_TADs_GMAP/"
outputFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/27_Number_TADs_per_Chromosome/table.txt"

# Execution
inFile.write(" ".join([prefix, categoryList, inputLocation, outputFileName])+"\n")

inFile.close()


