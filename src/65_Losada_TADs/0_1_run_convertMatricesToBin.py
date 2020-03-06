
# Import
import os
import sys

# Matrix List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/36_Losada_TADs/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/1_Sparse_Matrices/25K_norm/"
matrixList = ["LCONT", "SISA1", "SISA2"]

# Opening file
inFileName = fl + "01_cmb.txt"
inFile = open(inFileName,"w")

# Matrix Loop
for matrix in matrixList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Input
    chromosome = chrom
    resolution = "25000"
    matrixFileName = il + matrix + "/" + chrom + ".txt"
    chromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
    outputFileName = il + matrix + "/" + chrom + "_bin.txt"

    # Execution
    #command = "python 0_convertMatricesToBin.py "+" ".join([chromosome, resolution, matrixFileName, chromSizesFileName, outputFileName])
    #os.system(command)
    inFile.write(" ".join([chromosome, resolution, matrixFileName, chromSizesFileName, outputFileName])+"\n")

inFile.close()


