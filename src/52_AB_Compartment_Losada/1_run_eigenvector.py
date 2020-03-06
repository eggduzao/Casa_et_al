
# Import
import os
import sys

# Condition List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/23_AB_Compartment_Losada/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/19_Process_Losada_Hic/2_Merged_Hic/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/23_AB_Compartment_Losada/1_Eigenvector/"
condList = ["LCONT", "SISA1", "SISA2"]
condNameList = ["LCONT", "SISA1", "SISA2"]

# Opening file
inFileName = fl + "1_eig.txt"
inFile = open(inFileName,"w")

# Condition Loop
for i in range(0,len(condList)):

  cond = condList[i]
  condName = condNameList[i]

  # Chromosome List
  chromList = [str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Name
    outLoc = ol + condName + "/"
    os.system("mkdir -p "+outLoc)
    name = "chr" + chrom

    # Input
    normType = "KR"
    inputFileName = il + cond + "/inter_30.hic"
    chromN = chrom
    frag = "BP"
    resolution = "100000"
    outputFileName = outLoc + name + ".txt"

    # Write
    inFile.write(" ".join([normType, inputFileName, chromN, frag, resolution, outputFileName])+"\n")

inFile.close()


