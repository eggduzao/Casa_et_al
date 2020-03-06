
# Import
import os
import sys
from glob import glob

# Bam List
counter = 1
chromSizesFile = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/input/"
il = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/2_Merged_Bam_Files/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/5_BigWig/"
bamFileNameList = ["SA1_MCF10A_CHIPSEQ", "SA2_MCF10A_CHIPSEQ", "SMC1_MCF10A_CHIPSEQ", "ZMYM2_MCF10A_CHIPSEQ", "Input_MCF10A_CHIPSEQ", "INPUT_MCF10A_CHIPSEQ_siSA2", "INPUT_MCF10A_CHIPSEQ_siSA1", "INPUT_MCF10A_CHIPSEQ_Control", "SA1_MCF10A_CHIPSEQ_Control", "SA2_MCF10A_CHIPSEQ_Control", "SA1_MCF10A_CHIPSEQ_siSA1", "SA2_MCF10A_CHIPSEQ_siSA1", "SA1_MCF10A_CHIPSEQ_siSA2", "SA2_MCF10A_CHIPSEQ_siSA2", "CTCF_MCF10A_CHIPSEQ", "INPUT_CTCF_MCF10A_CHIPSEQ"]

# Opening file
inputFileName = fl + "6_cbw.txt"
inputFile = open(inputFileName, "w")

# Bam Loop
for bamName in bamFileNameList:

  # Input
  downstreamExt = "0"
  upstreamExt = "200"
  genomeSizesFileName = chromSizesFile
  bamFileName = il + bamName + ".bam"
  bwFileName = ol + bamName + ".bw"

  # Execution
  inputFile.write(" ".join([downstreamExt, upstreamExt, genomeSizesFileName, bamFileName, bwFileName])+"\n")
  counter += 1

# Close
inputFile.close()


