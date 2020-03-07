
# Import
import os
import sys

# Hic List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/39_Process_Rao_Matrix/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/39_Process_Rao_Matrix/"
tl = "/usr/users/egadegu/scratch/SPL/"
hicList = ["Untreated_5000", "Untreated_10000", "Untreated_25000", "Untreated_50000", "Untreated_100000", "Untreated_250000", "Untreated_500000", 
           "Untreated_Synchronized_5000", "Untreated_Synchronized_10000", "Untreated_Synchronized_25000", "Untreated_Synchronized_50000", "Untreated_Synchronized_100000", "Untreated_Synchronized_250000", "Untreated_Synchronized_500000"]

# Opening file
inFileName = fl + "2_spl.txt"
inFile = open(inFileName,"w")

# Fasta Loop
for hicName in hicList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in list(range(1,23))+["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Input
    chromosome = chrom
    inputFileName = il + hicName + ".txt"
    temporaryLocation = tl + hicName + "/" + chrom + "/"
    outputFileName = il + hicName + "/" + chrom + ".txt"

    # Creating files
    inFile.write(" ".join([chromosome, inputFileName, temporaryLocation, outputFileName])+"\n")

# Closing file
inFile.close()

