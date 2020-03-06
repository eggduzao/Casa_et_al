
# Import
import os
import sys

# Folder list
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/1_Sparse_Matrices/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/4_TAD_Interaction_Decay/"
chromSFN = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
folderList = ["25K_none", "25K_norm", "50K_none", "50K_norm", "100K_none", "100K_norm", "250K_none", "250K_norm"]

# Open file
inputFileName = fl + "4_tid.txt"
inFile = open(inputFileName, "w")

# Folder loop
for folder in folderList:

  # Matrix list
  matrixPairList = [["SISA1", "LCONT"], ["SISA2", "LCONT"]]

  # Matrix loop
  for matrixName in matrixPairList:

    # Chromosome list
    chromList = ["chr" + str(e) for e in range(1,23)+["X"]]

    # Chromosome loop
    for chrom in chromList:

      # Name
      mName = matrixName[0]
      name = folder+"/"+mName+"/"+chrom
      nameM = folder+"/"+matrixName[0]+"/"+chrom
      nameP = folder+"/"+matrixName[1]+"/"+chrom

      # Parameters
      chromosome = chrom
      if("25K" in folder): resolution = "25000"
      elif("50K" in folder): resolution = "50000"
      elif("100K" in folder): resolution = "100000"
      elif("250K" in folder): resolution = "250000"
      chromSizesFileName = chromSFN
      matrixMinusFileName = il+nameM+".txt"
      matrixPlusFileName = il+nameP+".txt"
      outputFileName = ol+name+".txt"

      # Execution
      inFile.write(" ".join([chromosome, resolution, chromSizesFileName, matrixMinusFileName, matrixPlusFileName, outputFileName])+"\n")

inFile.close()
 

