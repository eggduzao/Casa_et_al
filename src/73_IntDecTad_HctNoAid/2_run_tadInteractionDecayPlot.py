
# Import
import os
import sys

# Folder list
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/44_FigS9C_IntDecTad_RaoNeg/1_TAD_Interaction_Decay/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/44_FigS9C_IntDecTad_RaoNeg/2_TAD_Interaction_Decay_Plot/"
folderList1 = ["STAG1_v_Untreated"]
folderList2 = ["STAG2_v_Untreated"]

# Folder loop
for i in range(0,len(folderList1)):

  # Chromosome list
  chromList = ["chr"+str(e) for e in list(range(1,23))+["X"]]
  #           1      2       3        4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21      22      X
  msa1List = ["1.0", "-0.5", "-0.98", "-0.5", "-1.2", "-0.2", "-0.4", "-0.6", "-0.7", "-1.0", "-1.0", "-1.0", "-1.3", "-1.4", "-0.8", "-0.8", "-0.3", "-1.0", "-1.0", "-0.8", "-0.5", "-1.0", "-0.8"]
  msa2List = ["1.0", "-0.5", "-0.77", "-0.5", "-1.0", "+0.0", "-0.2", "-0.6", "-0.7", "-1.0", "-0.8", "-0.9", "-1.2", "-1.4", "-0.8", "-0.8", "-0.2", "-1.0", "-1.0", "-0.8", "-0.5", "-1.0", "-0.8"]

  # Chromosome loop
  for j in range(0,len(chromList)):

    # Name
    chrom = chromList[j]

    # Input
    chromosome = chrom
    resolution = "50000"
    msa1 = msa1List[j]
    msa2 = msa2List[j]
    inputSA1TableFileName = il + folderList1[i] + "/" + chrom + ".txt"
    inputSA2TableFileName = il + folderList2[i] + "/" + chrom + ".txt"
    outputFileName = ol + chrom + ".pdf"

    # Execution
    command = "R CMD BATCH '--args '"+chromosome+"' '"+resolution+"' '"+msa1+"' '"+msa2+"' '"+inputSA1TableFileName+"' '"+inputSA2TableFileName+"' '"+outputFileName+" /usr/users/egadegu/Projects/Wendt_Stag/Code/44_FigS9C_IntDecTad_RaoNeg/2_tadInteractionDecayPlot.R /usr/users/egadegu/Projects/Wendt_Stag/Code/44_FigS9C_IntDecTad_RaoNeg/2_tadInteractionDecayPlot.Rout"
    os.system(command)


