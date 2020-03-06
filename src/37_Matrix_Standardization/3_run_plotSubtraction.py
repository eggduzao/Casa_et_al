
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Matrix List
chromSizesFile = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
ml = "/home/egg/Projects/Wendt_Stag/Results/6_Matrix_Standardization/2_subtractions/250K/"
olm = "/home/egg/Projects/Wendt_Stag/Results/6_Matrix_Standardization/3_plots/matrix/"
olp = "/home/egg/Projects/Wendt_Stag/Results/6_Matrix_Standardization/3_plots/plot/"
matrixList = ["STAG1_plus_auxin_binarized__MINUS__STAG1_minus_auxin_binarized", "STAG2_plus_auxin_binarized__MINUS__STAG2_minus_auxin_binarized"]

# Matrix Loop
for matrix in matrixList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23) + ["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Parameters
    if("binarized" in matrix): ext = "bin"
    else: ext = "std"

    # Input
    chrom = chrom
    resolution = "250000"
    chromSizesFileName = chromSizesFile
    matrixFileName = ml + matrix + ".txt"
    outputMatrixFileName = olm + matrix.split("_")[0] + "_" + ext + "_" + chrom + ".txt"

    # Creating matrix
    command = "python 3_createFullTable.py "+" ".join([chrom, resolution, chromSizesFileName, matrixFileName, outputMatrixFileName])
    os.system(command)

    # Input
    minV = "-1"
    maxV = "1"
    outputPlotFileName = olp + matrix.split("_")[0] + "_"+ ext + "_" + chrom + ".pdf"

    # Creating plot
    command = "R CMD BATCH '--args '"+minV+"' '"+maxV+"' '"+outputMatrixFileName+"' '"+outputPlotFileName+" 3_hic.R 3_hic.Rout"
    os.system(command)


