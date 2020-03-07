
# Import
import os
import sys

# Folder list
chromSizesFile = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/44_FigS9C_IntDecTad_RaoNeg/input/"
ilt = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/25K_norm/"
ilc = "/usr/users/egadegu/Projects/Wendt_Stag/Results/39_Process_Rao_Matrix/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/44_FigS9C_IntDecTad_RaoNeg/1_TAD_Interaction_Decay/"

# Open file
inputFileName = fl + "1_tid.txt"
inFile = open(inputFileName, "w")

# Matrix list
matrixPairList = [[ilt+"69_127plus", ilc+"Untreated_25000"], [ilt+"3B9_5plus", ilc+"Untreated_25000"], [ilt+"69_127plus", ilc+"Untreated_Synchronized_25000"], [ilt+"3B9_5plus", ilc+"Untreated_Synchronized_25000"]]

# Matrix loop
for matrixPair in matrixPairList:

  # Chromosome list
  chromList = ["chr" + str(e) for e in list(range(1,23))+["X"]]

  # Chromosome loop
  for chrom in chromList:

    # Name
    tName = matrixPair[0].split("/")[-1]
    cName = matrixPair[1].split("/")[-1]
    if(tName == "69_127plus"): tName = "STAG1+"
    else: tName = "STAG2+"
    if(cName == "Untreated_25000"): cName = "Untreated"
    else: cName = "Synchronized"

    # Parameters
    chromosome = chrom
    resolution = "25000"
    chromSizesFileName = chromSizesFile
    matrixMinusFileName = matrixPair[1] + "/" + chrom + ".txt"
    matrixPlusFileName = matrixPair[0] + "/"+ chrom + ".txt"
    outputFileName = ol + "_".join([tName, "v", cName]) + "/" + chrom + ".txt"

    # Execution
    inFile.write(" ".join([chromosome, resolution, chromSizesFileName, matrixMinusFileName, matrixPlusFileName, outputFileName])+"\n")

inFile.close()
 

