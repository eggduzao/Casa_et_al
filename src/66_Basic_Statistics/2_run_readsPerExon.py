
# Import
import os
import sys

# Bam List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/37_Basic_Statistics/input/"
pfn = "/usr/users/egadegu/Projects/Wendt_Stag/Results/2_Genomic_Distribution/1_genomic_regions/regions.bed"
bl = "/usr/users/egadegu/Projects/Wendt_Stag/Data/expression/bam/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/37_Basic_Statistics/2_RNAseq/"
bamList = ["STAG1_aux_rep1", "STAG1_aux_rep2", "STAG1_noaux_rep1", "STAG1_noaux_rep2", "STAG2_aux_rep1", "STAG2_aux_rep2", "STAG2_noaux_rep1", "STAG2_noaux_rep2"]

# Open File
inputFileName = fl + "2_rpe.txt"
inFile = open(inputFileName, "w")

# Bam Loop
for bamName in bamList:

  # Input
  regionFileName = pfn
  bamFileName = bl + bamName + ".bam"
  outputFileName = ol + bamName + ".txt"

  # Execution
  inFile.write(" ".join([regionFileName, bamFileName, outputFileName])+"\n")

inFile.close()


