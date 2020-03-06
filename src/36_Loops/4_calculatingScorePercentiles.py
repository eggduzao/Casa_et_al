
# Import
import os
import sys
from numpy import isnan, isfinite, array, percentile

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
  inputFileName = inputLocation + chrom + ".txt.gz"
  uncompressedFileName = tempLoc + "uncompressedFileName.txt"
  command = "gzip -cd "+inputFileName+" > "+uncompressedFileName
  os.system(command)

  # Updating vector
  inputContactFile = open(uncompressedFileName, "rU")
  inputContactFile.readline()
  for line in inputContactFile:
    ll = line.strip().split("\t")
    score = float(ll[5])
    if(isnan(score) or (not isfinite(score))): continue
    scoreVec.append(score)
  inputContactFile.close()

# Fetching percentiles
scoreVec = array(scoreVec)
outputFile = open(outputFileName, "w")
for i in percList: outputFile.write("\t".join([str(i), str(percentile(scoreVec, i))])+"\n")
outputFile.close()

# Removing temp
command = "rm -rf "+tempLoc
os.system(command)


