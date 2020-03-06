
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import itertools
import numpy as np
from functools import reduce

# Input
resolution = int(sys.argv[1])
allowedShiftBins = int(sys.argv[2])
tadLabelList = sys.argv[3].split(",")
tadFileNameList = sys.argv[4].split(",")
tempLoc = sys.argv[5]
outputFileName = sys.argv[6]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

# Keeping only first order TADs
def keepFirstOnly(inFileName, outFileName):

  inFile = open(inFileName, "rU")
  outFile = open(outFileName, "w")
  for line in inFile:
    ll = line.strip().split("\t")
    if(ll[4] != "1"): continue
    outFile.write("\t".join(ll[:3])+"\n")
  inFile.close()
  outFile.close()

# Fetching list of possible TADs
def fileToList(inFileName):

  vector = []
  inFile = open(inFileName, "rU")
  for line in inFile:
    ll = line.strip().split("\t")
    start = int(ll[1])
    end = int(ll[2])
    for i in range(start-(resolution*allowedShiftBins), start+(resolution*allowedShiftBins)+1, resolution):
      for j in range(end-(resolution*allowedShiftBins), end+(resolution*allowedShiftBins)+1, resolution):
        vector.append(":".join([ll[0], str(i), str(j)]))
  inFile.close()
  return vector

# Overlap function
def overlap(matrix):
  result = reduce(np.intersect1d, matrix)
  return len(result)

###################################################################################################
# Execution
###################################################################################################

# Keeping only first order TADs
foTadFileNameList = []
for i in range(0,len(tadLabelList)):
  tadFileName = tadFileNameList[i]
  label = tadLabelList[i]
  tempFileName = tempLoc + label + "_fo.bed"
  keepFirstOnly(tadFileName, tempFileName)
  foTadFileNameList.append(tempFileName)

print foTadFileNameList

# Fetching vectors of each file
vectorList = []
for i in range(0,len(foTadFileNameList)):
  foTadFileName = foTadFileNameList[i]
  vector = fileToList(foTadFileName)
  vectorList.append(vector)

# Writing overlap
l = len(tadLabelList)
combList = list(itertools.combinations(range(l), 1)) + list(itertools.combinations(range(l), 2)) + list(itertools.combinations(range(l), 3)) + list(itertools.combinations(range(l), 4))
outputFile = open(outputFileName, "w")
for combL in combList:
  inLabels = []
  inMatrix = []
  for k in combL:
    inLabels.append(tadLabelList[k])
    inMatrix.append(vectorList[k])
  ovl = overlap(inMatrix)
  outputFile.write("\t".join([" + ".join(inLabels), str(ovl)])+"\n")
outputFile.close()


