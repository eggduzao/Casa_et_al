
# Import
import os
import sys

# Folder list
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/43_FigS9B_IntDec_RaoNeg/1_Interaction_Decay/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/43_FigS9B_IntDec_RaoNeg/2_Interaction_Decay_Plot/"
treatList = ["STAG1+", "STAG2+"]
contrList = ["Synchronized", "Synchronized"]

# Folder loop
for i in range(0,len(treatList)):

  # Chromosome list
  treatName = treatList[i]
  contrName = contrList[i]
  chromList = ["chr"+str(e) for e in list(range(1,23))+["X"]]
  if(treatName == "STAG1+"): ss = "0.4"
  else: ss = "0.28"

  # Chromosome loop
  for chrom in chromList:

    # Name
    outputLocation = ol + "_v_".join([treatName, contrName]) +"/"
    os.system("mkdir -p "+outputLocation)

    # Parameters
    binSize = "50"
    suby = ss
    chromosome = chrom
    cName = "Untreated"
    tName = treatName
    inputTableFileName1 = il + contrName + "/" + chrom + ".txt"
    inputTableFileName2 = il + treatName + "/" + chrom + ".txt"
    outputFileName = outputLocation + chrom + ".pdf"

    # Execution
    command = "R CMD BATCH '--args '"+binSize+"' '"+suby+"' '"+chromosome+"' '"+cName+"' '"+tName+"' '"+inputTableFileName1+"' '"+inputTableFileName2+"' '"+outputFileName+" /usr/users/egadegu/Projects/Wendt_Stag/Code/43_FigS9B_IntDec_RaoNeg/2_interactionDecayPlot.R /usr/users/egadegu/Projects/Wendt_Stag/Code/43_FigS9B_IntDec_RaoNeg/2_interactionDecayPlot.Rout"
    os.system(command)


