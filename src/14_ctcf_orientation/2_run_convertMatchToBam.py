
# Import
import os
import sys
from glob import glob

# Bed file list
counter = 1
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/14_ctcf_orientation/input_b2b/"
bl = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/tf/motifs/"
bedList = glob(bl+"*.bed")

# Bed Loop
for bedFile in bedList:

  # Parameters
  name = bedFile.split("/")[-1].split(".")[0]

  # Input 
  bedFileName = bedFile
  chromSizesFileName = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
  tempLoc = "/scratch/eduardo/"+name+"/"
  outputFileName = bl+name+".bam"

  # Input File
  inFileName = il+str(counter)+".txt"
  inFile = open(inFileName, "w")
  inFile.write("\n".join([bedFileName, chromSizesFileName, tempLoc, outputFileName]))
  inFile.close()
  counter += 1


