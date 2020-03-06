
# Import
import os
import sys

# Bam List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/18_Process_Losada_ChIP/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/2_Merged_Bam_Files/"
tl = "/scratch/egadegu/MAC/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/4_Peaks/"

treatList = ["SMC1_MCF10A_ChIP-seq", "ZMYM2_MCF10A_ChIP-seq", "SA1_MCF10A_ChIP-seq", "SA2_MCF10A_ChIP-seq",
             "SA1_MCF10A_ChIP-seq_Control", "SA2_MCF10A_ChIP-seq_Control", "SA1_MCF10A_ChIP-seq_siSA1", "SA1_MCF10A_ChIP-seq_siSA2", "SA2_MCF10A_ChIP-seq_siSA1", "SA2_MCF10A_ChIP-seq_siSA2"]
controlList = ["Input_MCF10A_ChIP-seq", "Input_MCF10A_ChIP-seq", "Input_MCF10A_ChIP-seq", "Input_MCF10A_ChIP-seq",
               "INPUT_MCF10A_ChIP-seq_Control", "INPUT_MCF10A_ChIP-seq_Control", "INPUT_MCF10A_ChIP-seq_siSA1", "INPUT_MCF10A_ChIP-seq_siSA2", "INPUT_MCF10A_ChIP-seq_siSA1", "INPUT_MCF10A_ChIP-seq_siSA2"]

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


