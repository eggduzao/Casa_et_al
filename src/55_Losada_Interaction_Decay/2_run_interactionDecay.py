
# Import
import os
import sys

# Folder list
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/1_Sparse_Matrices/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/2_Interaction_Decay/"
folderList = ["25K_none", "25K_norm", "50K_none", "50K_norm", "100K_none", "100K_norm", "250K_none", "250K_norm"]

# Open file
inputFileName = fl + "2_rid.txt"
inFile = open(inputFileName, "w")

# Folder loop
for folder in folderList:

  # Matrix list
  matrixList = ["LCONT", "SISA1", "SISA2"]

  # Matrix loop
  for matrixName in matrixList:

    # Chromosome list
    chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

    # Chromosome loop
    for chrom in chromList:

      # Name
      name = folder + "/" + matrixName + "/" + chrom
      #print name

      # Parameters
      chromosome = chrom
      meanHalfSize = "100"
      if("25K" in folder): resolution = "25000"
      elif("50K" in folder): resolution = "50000"
      elif("100K" in folder): resolution = "100000"
      elif("250K" in folder): resolution = "250000"
      chromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
      matrixFileName = il + name + ".txt"
      outputFileName = ol + name + ".txt"

      # Execution
      #command = "python /usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/2_interactionDecay.py "+" ".join([chromosome, meanHalfSize, resolution, chromSizesFileName, matrixFileName, outputFileName])
      #os.system(command)
      inFile.write(" ".join([chromosome, meanHalfSize, resolution, chromSizesFileName, matrixFileName, outputFileName])+"\n")

inFile.close()


