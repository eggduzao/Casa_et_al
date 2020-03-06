
# Import
import os
import sys

# Bam List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/18_Process_Losada_ChIP/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/1_Raw_Bam_Files/"
tl = "/scratch/egadegu/MEB/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/2_Merged_Bam_Files/"
bamListList = [["Input_MCF10A_ChIP-seq"], ["INPUT_MCF10A_ChIP-seq_Control"], ["INPUT_MCF10A_ChIP-seq_siSA1"], ["INPUT_MCF10A_ChIP-seq_siSA2"], ["SMC1_MCF10A_ChIP-seq"], ["ZMYM2_MCF10A_ChIP-seq"], ["SA1_MCF10A_ChIP-seq"], ["SA2_MCF10A_ChIP-seq"], 
               ["SA1_MCF10A_ChIP-seq_Control_Rep_1", "SA1_MCF10A_ChIP-seq_Control_Rep_2"], ["SA1_MCF10A_ChIP-seq_siSA1_Rep_1", "SA1_MCF10A_ChIP-seq_siSA1_Rep_2"], ["SA1_MCF10A_ChIP-seq_siSA2_Rep_1", "SA1_MCF10A_ChIP-seq_siSA2_Rep_2"], 
               ["SA2_MCF10A_ChIP-seq_Control_Rep_1", "SA2_MCF10A_ChIP-seq_Control_Rep_2"], ["SA2_MCF10A_ChIP-seq_siSA1_Rep_1", "SA2_MCF10A_ChIP-seq_siSA1_Rep_2"], ["SA2_MCF10A_ChIP-seq_siSA2_Rep_1", "SA2_MCF10A_ChIP-seq_siSA2_Rep_2"]]

# Opening file
inFileName = fl + "2_meb.txt"
inFile = open(inFileName,"w")

# Bam Loop
for bamList in bamListList:

  # Parameters
  if("_Rep_1" in bamList[0]): bamName = bamList[0].split("_Rep_")[0]
  else: bamName = bamList[0]

  # Input
  bamFileNameList = ",".join([il + e + ".bam" for e in bamList])
  tempBamFileName = tl + bamName + "/" + bamName + ".bam"
  outputBamFileName = ol + bamName + ".bam"

  # Creating files
  inFile.write(" ".join([bamFileNameList, tempBamFileName, outputBamFileName])+"\n")

# Closing file
inFile.close()


