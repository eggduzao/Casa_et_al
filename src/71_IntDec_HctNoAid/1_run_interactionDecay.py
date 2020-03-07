
# Import
import os
import sys

# Matrix list
chromSizesFile = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/43_FigS9B_IntDec_RaoNeg/input/"
il1 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/50K_norm/"
il2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/39_Process_Rao_Matrix/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/43_FigS9B_IntDec_RaoNeg/1_Interaction_Decay/"
#matrixFileNameList = [il1+"3B9_5-", il1+"69_127-", il2+"Untreated_50000", il2+"Untreated_Synchronized_50000"]
#outNameList = ["STAG2-", "STAG1-", "Untreated", "Synchronized"]
matrixFileNameList = [il1+"3B9_5plus", il1+"69_127plus"]
outNameList = ["STAG2+", "STAG1+"]

# Open file
inputFileName = fl + "1_cid.txt"
inFile = open(inputFileName, "w")

# Matrix loop
for i in range(0,len(matrixFileNameList)):

  # Chromosome list
  matrixFName = matrixFileNameList[i]
  outName = outNameList[i]
  chromList = ["chr"+str(e) for e in list(range(1,23))+["X"]]

  # Chromosome loop
  for chrom in chromList:

    # Parameters
    chromosome = chrom
    meanHalfSize = "100"
    resolution = "50000"
    chromSizesFileName = chromSizesFile
    matrixFileName = matrixFName + ".txt"
    outputFileName = ol + outName + "/" + chrom + ".txt"

    # Execution
    inFile.write(" ".join([chromosome, meanHalfSize, resolution, chromSizesFileName, matrixFileName, outputFileName])+"\n")

inFile.close()


