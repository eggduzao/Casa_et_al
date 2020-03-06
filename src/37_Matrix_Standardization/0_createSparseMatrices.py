
# Import
import os
import sys
import math
import numpy as np

# Input
juicerCommand = sys.argv[1]
kindOfMatrix = sys.argv[2]
kindOfNormalization = sys.argv[3]
unitOfResolution = sys.argv[4]
resolution = sys.argv[5]
chromSizesFileName = sys.argv[6]
inputHicFileName = sys.argv[7]
tempLocation = sys.argv[8]
outputMatrixFileName = sys.argv[9]
outputScoreFileName = sys.argv[10]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)
hicName = inputHicFileName.split("/")[-1].split(".")[0]

# Chromosome sizes dict
chromSizeDict = dict()
chromSizesFile = open(chromSizesFileName,"rU")
for line in chromSizesFile:
  ll = line.strip().split("\t")
  chromSizeDict[ll[0]] = ll[1]
chromSizesFile.close()

# Opening output file to merge everything on the go
outputMatrixFile = open(outputMatrixFileName,"w")
scoreVector = []

# Creating individual sparse chromosome matrices
for chrom in sorted(chromSizeDict.keys()):

  # Initialization
  chrWoChr = chrom.split("chr")[-1]
  region = ":".join([chrWoChr,"1",str(chromSizeDict[chrom])])

  # Creating sparse matrix
  tempFileName = tempLocation+chrom+".txt"
  command = " ".join([juicerCommand, "dump", kindOfMatrix, kindOfNormalization, inputHicFileName,
                      region, region, unitOfResolution, resolution, tempFileName])
  os.system(command)

  # Writing entries
  tempFile = open(tempFileName,"rU")
  for line in tempFile:
    ll = line.strip().split("\t")
    value = float(ll[2])
    if(math.isnan(value) or not np.isfinite(value)): continue
    scoreVector.append(value)
    outputMatrixFile.write("\t".join([chrom]+ll)+"\n")
  tempFile.close()

  # Tracking
  print "Step 1: ", hicName, chrom

# Calculating percentiles
outputScoreFile = open(outputScoreFileName,"w")
scoreVectorNp = np.array(scoreVector)
for i in range(0,101):
  outputScoreFile.write("\t".join([str(i),str(np.percentile(scoreVectorNp, i))])+"\n")
  print "Step 2: ", i
outputScoreFile.close()

# Termination
outputMatrixFile.close()
command = "rm -rf "+tempLocation
os.system(command)


