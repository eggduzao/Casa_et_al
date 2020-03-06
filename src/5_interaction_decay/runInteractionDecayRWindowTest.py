
# Import
import os
import sys

# Folder list
il = "/home/egg/Projects/Papantonis_Stag/Results/interaction_decay/tables/"
ol = "/home/egg/Projects/Papantonis_Stag/Results/interaction_decay/graphs/window_test/"
folderList = ["10K_norm", "25K_norm", "50K_norm"]

# Folder loop
for folder in folderList:

  # Matrix list
  matrixList = ["3B9_5"]

  # Matrix loop
  for matrixName in matrixList:

    # Chromosome list
    chromList = ["chr22"]

    # Chromosome loop
    for chrom in chromList:

      # Window list
      windowList = ["3", "5", "10", "20", "30", "40", "50", "75", "100", "125", "150", "175", "200"]

      # Window loop
      for window in windowList:
    
        # Name
        name = folder+"_bin"+window
        outputLocation = ol
        os.system("mkdir -p "+outputLocation)

        # Parameters
        binSize = window
        chromosome = chrom
        inputTableFileName1 = il+folder+"/"+matrixName+"-/"+chrom+".txt"
        inputTableFileName2 = il+folder+"/"+matrixName+"plus/"+chrom+".txt"
        outputFileName = ol+name+".png"

        # Execution
        command = "R CMD BATCH '--args '"+binSize+"' '"+chromosome+"' '"+inputTableFileName1+"' '"+inputTableFileName2+"' '"+outputFileName+" interactionDecay.R interactionDecay.Rout"
        os.system(command)


