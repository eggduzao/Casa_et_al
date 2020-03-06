
# Import
import os
import sys

# Matrix List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_matrix_files/10K_norm/"
matrixList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

# Matrix Loop
for matrix in matrixList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Input
    chromosome = chrom
    resolution = "10000"
    matrixFileName = il+matrix+"/"+chrom+".txt"
    chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
    outputFileName = il+matrix+"/"+chrom+"_bin.txt"

    # Execution
    print matrix, chrom
    command = "python 0_convertMatricesToBin.py "+" ".join([chromosome, resolution, matrixFileName, chromSizesFileName, outputFileName])
    os.system(command)


