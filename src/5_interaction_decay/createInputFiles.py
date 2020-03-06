
# Import
import os
import sys

# Folder list
counter = 1
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_matrix_files/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/interaction_decay/tables/"
folderList = ["10K_none", "10K_norm", "50K_none", "50K_norm"]
outLocation = "./input/"

# Folder loop
for folder in folderList:

  # Matrix list
  matrixList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

  # Matrix loop
  for matrixName in matrixList:

    # Chromosome list
    chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

    # Chromosome loop
    for chrom in chromList:

      # Name
      name = folder+"/"+matrixName+"/"+chrom
      print name

      # Parameters
      chromosome = chrom
      meanHalfSize = "100"
      if("10K" in folder): resolution = "10000"
      else: resolution = "50000"
      chromSizesFileName = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/chrom.sizes.hg19"
      matrixFileName = il+name+".txt"
      outputFileName = ol+name+".txt"

      # Execution
      inFile = open(outLocation+str(counter),"w")
      inFile.write("\n".join([chromosome, meanHalfSize, resolution, chromSizesFileName, matrixFileName, outputFileName]))
      inFile.close()
      counter += 1


