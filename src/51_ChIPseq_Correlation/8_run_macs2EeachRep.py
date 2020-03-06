
# Import
import os
import sys

# Bam List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/1_Raw_Bam_Files/"
tl = "/scratch/egadegu/MACR/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/8_Peaks_Separated_Replicate/"
treatList = ["SA1_Rep1", "SA1_Rep2", "SA2_Rep1", "SA2_Rep2"]
controlList = ["SA1_Input", "SA1_Input", "SA2_Input", "SA2_Input"]
threshList = ["0.005", "0.005", "0.006", "0.006"]

# Opening file
inFileName = fl + "8_macr.txt"
inFile = open(inFileName, "w")

# Bam Loop
for i in range(0,len(treatList)):

  # Parameters
  treatName = treatList[i]
  controlName = controlList[i]
  thresh = threshList[i]

  # Input
  threshold = thresh
  dataType = "CHIP"
  treatmentFileName = il + treatName + ".bam"
  controlFileName = il + controlName + ".bam"
  tempLocation = tl + treatName + "/"
  outputLocation = ol

  # Creating files
  inFile.write(" ".join([threshold, dataType, treatmentFileName, controlFileName, tempLocation, outputLocation])+"\n")

# Close
inFile.close()


