
# Import
import os
import sys

# Folder list
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/4_TAD_Interaction_Decay/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/5_TAD_Interaction_Decay_Plot/"
folderList = ["25K_none", "25K_norm", "50K_none", "50K_norm", "100K_none", "100K_norm", "250K_none", "250K_norm"]

# Folder loop
for folder in folderList:

  # Chromosome list
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome loop
  for chrom in chromList:

    # Name
    print folder, chrom
    outputLocation = ol+folder+"/"
    os.system("mkdir -p "+outputLocation)

    # Input
    chromosome = chrom
    if("25K" in folder): resolution = "25000"
    elif("50K" in folder): resolution = "50000"
    elif("100K" in folder): resolution = "100000"
    elif("250K" in folder): resolution = "250000"
    inputSA1TableFileName = il+folder+"/SISA1/"+chrom+".txt"
    inputSA2TableFileName = il+folder+"/SISA2/"+chrom+".txt"
    outputFileName = outputLocation+chrom+".pdf"

    # Execution
    command = "R CMD BATCH '--args '"+chromosome+"' '"+resolution+"' '"+inputSA1TableFileName+"' '"+inputSA2TableFileName+"' '"+outputFileName+" /usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/5_tadInteractionDecayPlot.R /usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/5_tadInteractionDecayPlot.Rout"
    os.system(command)


