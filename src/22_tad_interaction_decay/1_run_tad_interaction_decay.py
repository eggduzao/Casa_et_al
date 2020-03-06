
# Import
import os
import sys

# Folder list
counter = 1
fn = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/22_tad_interaction_decay/1_input/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_matrix_files/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/22_tad_interaction_decay/1_tad_interaction_decay/"
chromSFN = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
folderList = ["10K_norm", "25K_norm", "50K_norm"]

# Folder loop
for folder in folderList:

  # Matrix list
  matrixPairList = [["3B9_5-", "3B9_5plus"], ["69_127-", "69_127plus"]]

  # Matrix loop
  for matrixName in matrixPairList:

    # Chromosome list
    chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

    # Chromosome loop
    for chrom in chromList:

      # Name
      mName = matrixName[0][:-1]
      name = folder+"/"+mName+"/"+chrom
      nameM = folder+"/"+matrixName[0]+"/"+chrom
      nameP = folder+"/"+matrixName[1]+"/"+chrom

      # Parameters
      chromosome = chrom
      if("10K" in folder): resolution = "10000"
      elif("25K" in folder): resolution = "25000"
      elif("50K" in folder): resolution = "50000"
      chromSizesFileName = chromSFN
      matrixMinusFileName = il+nameM+".txt"
      matrixPlusFileName = il+nameP+".txt"
      outputFileName = ol+name+".txt"

      # Execution
      inFileName = fn+str(counter)+".txt"
      inFile = open(inFileName, "w")
      inFile.write("\n".join([chromosome, resolution, chromSizesFileName, matrixMinusFileName, matrixPlusFileName, outputFileName]))
      inFile.close()
      counter += 1

      # Execution
      #command = "python 1_tad_interaction_decay.py "+" ".join([chromosome, resolution, chromSizesFileName, matrixMinusFileName, matrixPlusFileName, outputFileName])
      #os.system(command) 
 

