
# Import
import os
import sys

# Folder list
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/2_Interaction_Decay/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/3_Interaction_Decay_Plot/"
folderList = ["25K_none", "25K_norm", "50K_none", "50K_norm", "100K_none", "100K_norm", "250K_none", "250K_norm"]

# Folder loop
for folder in folderList:

  # Matrix list
  matrixList = ["SISA1", "SISA2"]

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
      inputTableFileName1 = il+folder+"/LCONT/"+chrom+".txt"
      inputTableFileName2 = il+folder+"/"+matrixName+"/"+chrom+".txt"
      outputFileName = ol+name+".pdf"

      # Execution
      command = "R CMD BATCH '--args '"+binSize+"' '"+chromosome+"' '"+inputTableFileName1+"' '"+inputTableFileName2+"' '"+outputFileName+" /usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/3_interactionDecayPlot.R /usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/3_interactionDecayPlot.Rout"
      os.system(command)


