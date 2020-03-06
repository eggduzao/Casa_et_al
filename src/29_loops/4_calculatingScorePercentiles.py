
# Import
import os
import sys
import numpy as np

# Input
inputLocation = sys.argv[1]
tempLoc = sys.argv[2]
outputFileName = sys.argv[3]

# Initialization
percList = range(0,101)
command = "mkdir -p "+tempLoc
os.system(command)

# Chromosome List
chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Chromosome Loop
scoreVec = []
for chrom in chromList:

  # Creating temporary uncompressed file
  inputFileName = inputLocation + chrom + ".spline_pass1.res25000.significances.txt.gz"
  uncompressedFileName = tempLoc + "uncompressedFileName.txt"
  command = "gzip -cd "+inputFileName+" > "+uncompressedFileName
  os.system(command)

  # Updating vector
  inputContactFile = open(uncompressedFileName, "rU")
  inputContactFile.readline()
  for line in inputContactFile:
    ll = line.strip().split("\t")
    scoreVec.append(float(ll[6]))
  inputContactFile.close()

# Fetching percentiles
scoreVec = np.array(scoreVec)
outputFile = open(outputFileName, "w")
for i in percList: outputFile.write("\t".join([str(i), str(np.percentile(scoreVec, i))])+"\n")
outputFile.close()

# Removing temp
command = "rm -rf "+tempLoc
os.system(command)


