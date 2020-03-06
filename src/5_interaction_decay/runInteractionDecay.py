
# Import
import os
import sys

# Folder list
il = "/media/egg/sbc/AG_Papantonis/Eduardo/Papantonis_Stag/Data/stag_matrix_files/"
ol = "/home/egg/Projects/Papantonis_Stag/Results/5_interaction_decay/tables/"
folderList = ["10K_none", "10K_norm", "25K_none", "25K_norm", "50K_none", "50K_norm"]

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
      elif("25K" in folder): resolution = "25000"
      elif("50K" in folder): resolution = "50000"
      chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19"
      matrixFileName = il+name+".txt"
      outputFileName = ol+name+".txt"

      # Execution
      command = "python interactionDecay.py "+" ".join([chromosome,meanHalfSize,resolution,chromSizesFileName,matrixFileName,outputFileName])
      os.system(command)


