
# Import
import os
import sys

# Bam List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/2_Merged_Bam_Files/"
tl = "/scratch/egadegu/MAC/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/4_Peaks/"
treatList = ["SA1", "SA2"]
controlList = ["SA1_Input", "SA2_Input"]

# Opening file
inFileName = fl + "4_mac.txt"
inFile = open(inFileName, "w")

# Bam Loop
for i in range(0,len(treatList)):

  # Parameters
  treatName = treatList[i]
  controlName = controlList[i]

  # Input
  dataType = "CHIP"
  treatmentFileName = il + treatName + ".bam"
  controlFileName = il + controlName + ".bam"
  tempLocation = tl + treatName + "/"
  outputLocation = ol

  # Creating files
  inFile.write(" ".join([dataType, treatmentFileName, controlFileName, tempLocation, outputLocation])+"\n")

# Close
inFile.close()


