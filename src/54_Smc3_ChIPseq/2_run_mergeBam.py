
# Import
import os
import sys

# Bam List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/1_Raw_Bam_Files/"
tl = "/scratch/egadegu/MEB/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/2_Merged_Bam_Files/"
bamListList = [["Smc3_Input"], ["Smc3_Rep1", "Smc3_Rep2"]]

# Opening file
inFileName = fl + "2_meb.txt"
inFile = open(inFileName,"w")

# Bam Loop
for bamList in bamListList:

  # Parameters
  if("_Rep1" in bamList[0]): bamName = bamList[0].split("_Rep")[0]
  else: bamName = bamList[0]

  # Input
  bamFileNameList = ",".join([il + e + ".bam" for e in bamList])
  tempBamFileName = tl + bamName + "/" + bamName + ".bam"
  outputBamFileName = ol + bamName + ".bam"

  # Creating files
  inFile.write(" ".join([bamFileNameList, tempBamFileName, outputBamFileName])+"\n")

# Closing file
inFile.close()


