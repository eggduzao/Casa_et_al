
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Bam List
il = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/2_Merged_Bam_Files/"
tl = "/scratch/egusmao/STAG_CTR/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/3_Read_Counts/"
bamFileNameList = ["SA1_MCF10A_CHIPSEQ", "SA2_MCF10A_CHIPSEQ", "SMC1_MCF10A_CHIPSEQ", "ZMYM2_MCF10A_CHIPSEQ", "CTCF_MCF10A_CHIPSEQ", "INPUT_CTCF_MCF10A_CHIPSEQ", "Input_MCF10A_CHIPSEQ", "INPUT_MCF10A_CHIPSEQ_siSA2", "INPUT_MCF10A_CHIPSEQ_siSA1", "INPUT_MCF10A_CHIPSEQ_Control", "SA1_MCF10A_CHIPSEQ_Control", "SA2_MCF10A_CHIPSEQ_Control", "SA1_MCF10A_CHIPSEQ_siSA1", "SA2_MCF10A_CHIPSEQ_siSA1", "SA1_MCF10A_CHIPSEQ_siSA2", "SA2_MCF10A_CHIPSEQ_siSA2"]

# Output file
outputFileName = ol + "number_of_valid_reads.txt"
outputFile = open(outputFileName,"w")

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


