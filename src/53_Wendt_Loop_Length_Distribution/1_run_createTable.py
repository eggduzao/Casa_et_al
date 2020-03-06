
# Import
import os
import sys

# Folder List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/24_Wendt_Loop_Length_Distribution/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/5_Loops/3_contact_files_filtered/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/24_Wendt_Loop_Length_Distribution/"
loopList = ["STAG1_minus_auxin_150Kb_5Mb" , "STAG1_plus_auxin_150Kb_5Mb" , "STAG2_minus_auxin_150Kb_5Mb" , "STAG2_plus_auxin_150Kb_5Mb"]
nameList = ["STAG1-AUX", "STAG1+AUX", "STAG2-AUX", "STAG2+AUX"]

# Opening file
inFileName = fl + "1_ctd.txt"
inFile = open(inFileName,"w")

# Input
inputNameList = ",".join(nameList)
inputFileNameList = ",".join([il + e + ".txt" for e in loopList])
outputFileName = ol + "table.txt"

# Creating files
inFile.write(" ".join([inputNameList, inputFileNameList, outputFileName])+"\n")

# Closing file
inFile.close()


