
# Import
import os
import sys
from glob import glob

# Fastq List
counter = 1
fl = "./input_fqd/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/rao/sra/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/rao/fasta/"
sraList = glob(il+"*.sra")

# Fastq Loop
for sraFile in sraList:

  # Name
  name = sraFile.split("/")[-1].split(".")[0]

  # Input
  sraFileName = sraFile
  outputLocation = ol + name + "/"

  # Creating files
  inFileName = fl+str(counter)+".txt"
  inFile = open(inFileName,"w")
  inFile.write("\n".join([sraFileName, outputLocation]))
  inFile.close()
  counter += 1


