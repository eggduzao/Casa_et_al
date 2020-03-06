
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Bam List
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/2_Merged_Bam_Files/"
tl = "/scratch/egadegu/CTR/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/3_Read_Counts/"
#bamFileNameList = ["Input_MCF10A_ChIP-seq", "INPUT_MCF10A_ChIP-seq_Control", "INPUT_MCF10A_ChIP-seq_siSA1", "INPUT_MCF10A_ChIP-seq_siSA2", "SMC1_MCF10A_ChIP-seq", "ZMYM2_MCF10A_ChIP-seq", "SA1_MCF10A_ChIP-seq", "SA2_MCF10A_ChIP-seq", "SA1_MCF10A_ChIP-seq_Control", "SA1_MCF10A_ChIP-seq_siSA1", "SA1_MCF10A_ChIP-seq_siSA2", "SA2_MCF10A_ChIP-seq_Control", "SA2_MCF10A_ChIP-seq_siSA1", "SA2_MCF10A_ChIP-seq_siSA2"]
bamFileNameList = ["SA1_MCF10A_ChIP-seq_Control", "SA1_MCF10A_ChIP-seq_siSA1", "SA1_MCF10A_ChIP-seq_siSA2", "SA2_MCF10A_ChIP-seq_Control", "SA2_MCF10A_ChIP-seq_siSA1", "SA2_MCF10A_ChIP-seq_siSA2"]

# Output file
command = "mkdir -p " + ol
os.system(command)
outputFileName = ol + "number_of_valid_reads2.txt"
outputFile = open(outputFileName, "w")

# Create temporary location
command = "mkdir -p "+tl
os.system(command)

# Bam Loop
for bamFileName in bamFileNameList:

  # Execution
  tempFileName = tl + "tempFileName.txt"
  os.system("samtools view -c "+ il + bamFileName + ".bam > " + tempFileName)
  tempFile = open(tempFileName,"rU")
  value = tempFile.readline().strip()
  tempFile.close()
  outputFile.write("\t".join([bamFileName,value])+"\n")
  os.system("rm "+tempFileName)

outputFile.close()


