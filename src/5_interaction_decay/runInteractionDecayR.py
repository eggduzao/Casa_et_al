
# Import
import os
import sys

# Folder list
il = "/home/egg/Projects/Papantonis_Stag/Results/5_interaction_decay/tables/"
ol = "/home/egg/Projects/Papantonis_Stag/Results/5_interaction_decay/graphs/"
folderList = ["10K_none", "10K_norm", "25K_none", "25K_norm", "50K_none", "50K_norm"]

# Folder loop
for folder in folderList:

  # Matrix list
  matrixList = ["3B9_5", "69_127"]

  # Matrix loop
  for matrixName in matrixList:

    # Chromosome list
    chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

    # Chromosome loop
    for chrom in chromList:

      # Name
      name = folder+"/"+matrixName+"/"+chrom
      outputLocation = ol+folder+"/"+matrixName+"/"
      os.system("mkdir -p "+outputLocation)

      # Parameters
      binSize = "50"
      chromosome = chrom
      inputTableFileName1 = il+folder+"/"+matrixName+"-/"+chrom+".txt"
      inputTableFileName2 = il+folder+"/"+matrixName+"plus/"+chrom+".txt"
      outputFileName = ol+name+".pdf"

      # Execution
      command = "R CMD BATCH '--args '"+binSize+"' '"+chromosome+"' '"+inputTableFileName1+"' '"+inputTableFileName2+"' '"+outputFileName+" interactionDecay.R interactionDecay.Rout"
      os.system(command)


