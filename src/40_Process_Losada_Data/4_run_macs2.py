
# Import
import os
import sys

# Bam List
counter = 1
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/input/"
il = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/2_Merged_Bam_Files/"
tl = "/scratch/egusmao/STAG_MAC/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/4_Peaks/"
treatList = ["SA1_MCF10A_CHIPSEQ", "SA2_MCF10A_CHIPSEQ", "SMC1_MCF10A_CHIPSEQ", "ZMYM2_MCF10A_CHIPSEQ", "SA1_MCF10A_CHIPSEQ_Control", "SA2_MCF10A_CHIPSEQ_Control", "SA1_MCF10A_CHIPSEQ_siSA1", "SA2_MCF10A_CHIPSEQ_siSA1", "SA1_MCF10A_CHIPSEQ_siSA2", "SA2_MCF10A_CHIPSEQ_siSA2", "CTCF_MCF10A_CHIPSEQ"]
controlList = ["Input_MCF10A_CHIPSEQ", "Input_MCF10A_CHIPSEQ", "Input_MCF10A_CHIPSEQ", "Input_MCF10A_CHIPSEQ", "INPUT_MCF10A_CHIPSEQ_Control", "INPUT_MCF10A_CHIPSEQ_Control", "INPUT_MCF10A_CHIPSEQ_siSA1", "INPUT_MCF10A_CHIPSEQ_siSA1", "INPUT_MCF10A_CHIPSEQ_siSA2", "INPUT_MCF10A_CHIPSEQ_siSA2", "INPUT_CTCF_MCF10A_CHIPSEQ"]

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
  tempLocation = tl + str(counter) + "/"
  outputLocation = ol

  # Creating files
  inFile.write(" ".join([dataType, treatmentFileName, controlFileName, tempLocation, outputLocation])+"\n")
  counter += 1

# Close
inFile.close()


