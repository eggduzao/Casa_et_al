
# Import
import os
import sys
from glob import glob

# Bam List
chromSizesFile = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.enhanced"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/2_Merged_Bam_Files/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/5_BigWig/"
bamFileNameList = ["Smc3_Input", "Smc3"]

# Opening file
inputFileName = fl + "6_cbw.txt"
inputFile = open(inputFileName, "w")

# Bam Loop
for bamName in bamFileNameList:

  # Input
  downstreamExt = "200"
  upstreamExt = "0"
  genomeSizesFileName = chromSizesFile
  bamFileName = il + bamName + ".bam"
  bwFileName = ol + bamName + ".bw"

  # Execution
  inputFile.write(" ".join([downstreamExt, upstreamExt, genomeSizesFileName, bamFileName, bwFileName])+"\n")

# Close
inputFile.close()


