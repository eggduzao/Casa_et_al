
# Import
import os
import sys

# Input List
chromSizesFile = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_matrix_files/25K_norm/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/fithic_matrix_files/25K_norm/"
os.system("mkdir -p "+ol)
inList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

# Input Loop
for inName in inList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Input
    resolution = "25000"
    chromSizesFileName = chromSizesFile
    inputMatrixFileName = il + inName + "/" + chrom + ".txt"
    outputFragmentFileName = ol + inName + "_" + chrom + "_fragment"
    outputContactFileName = ol + inName + "_" + chrom + "_contact"

    # Execution
    print inName, chrom
    command = "python 1_createInputFiles.py "+" ".join([resolution, chromSizesFileName, inputMatrixFileName, outputFragmentFileName, outputContactFileName])
    os.system(command)


