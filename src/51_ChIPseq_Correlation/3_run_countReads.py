
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Bam List
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/2_Merged_Bam_Files/"
tl = "/scratch/egadegu/CTR/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/3_Read_Counts/"
bamFileNameList = ["SA1_Input", "SA2_Input", "SA1", "SA2"]

# Output file
command = "mkdir -p " + ol
os.system(command)
outputFileName = ol + "number_of_valid_reads.txt"
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


