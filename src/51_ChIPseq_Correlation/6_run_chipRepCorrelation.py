
# Import
import os
import sys

# Bam List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/"
regionFile = "/usr/users/egadegu/rgtdata/hg19/binned_genome_100.bed"
wbl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/1_Raw_Bam_Files/"
lbl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/1_Raw_Bam_Files/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/6_Replicate_Correlation/"
treatList = [wbl+"SA1_Input.bam", wbl+"SA1_Rep1.bam", wbl+"SA1_Rep2.bam", wbl+"SA2_Input.bam", wbl+"SA2_Rep1.bam", wbl+"SA2_Rep2.bam", lbl+"Input_MCF10A_ChIP-seq.bam", lbl+"SA1_MCF10A_ChIP-seq.bam", lbl+"SA2_MCF10A_ChIP-seq.bam"]

# Opening file
inFileName = fl + "6_crc.txt"
inFile = open(inFileName, "w")

# Input
regionFileName = regionFile
bamFileNameList = ",".join(treatList)
outputFilePrefix = ol + "table"

# Creating files
inFile.write(" ".join([regionFileName, bamFileNameList, outputFilePrefix])+"\n")

# Close
inFile.close()


