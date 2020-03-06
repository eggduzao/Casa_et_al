
# Import
import os
import sys

# Folder list
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/22_tad_interaction_decay/1_tad_interaction_decay/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/22_tad_interaction_decay/2_tad_interaction_decay_plots/"
folderList = ["10K_norm", "25K_norm", "50K_norm"]
folderList = ["50K_norm"]

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
    if("10K" in folder): resolution = "10000"
    elif("25K" in folder): resolution = "25000"
    elif("50K" in folder): resolution = "50000"
    inputSA1TableFileName = il+folder+"/69_127/"+chrom+".txt"
    inputSA2TableFileName = il+folder+"/3B9_5/"+chrom+".txt"
    outputFileName = outputLocation+chrom+".pdf"

    # Execution
    command = "R CMD BATCH '--args '"+chromosome+"' '"+resolution+"' '"+inputSA1TableFileName+"' '"+inputSA2TableFileName+"' '"+outputFileName+" 1_tad_interaction_decay.R 1_tad_interaction_decay.Rout"
    os.system(command)


