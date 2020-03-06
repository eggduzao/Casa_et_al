
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import math

# Input
loopBins = int(sys.argv[1])
resolution = int(sys.argv[2])
regionFileName = sys.argv[3]
matrixFileName = sys.argv[4]
outputMatrixFileName = sys.argv[5]

###################################################################################################
# Functions
###################################################################################################

def readRegions(inputFileName):

  regions = []
  inputFile = open(inputFileName,"rU")
  for line in inputFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])
    region = [chrom, str(p1-(loopBins*resolution)), str(p1+((loopBins+1)*resolution)), chrom, str(p2-(loopBins*resolution)), str(p2+((loopBins+1)*resolution))]
    regions.append(region)
  inputFile.close()
  return regions

def readMatrix(matrixFileName):

  matrixDict = dict()
  matrixFile = open(matrixFileName,"rU")
  for line in matrixFile:
    ll = line.strip().split("\t")
    matrixDict[":".join(ll[:3])] = float(ll[3])
  matrixFile.close()

  return matrixDict

def updateMatrix(matrix, region1, region2, matrixDict, resolution):

  counterI = 0
  chrom = region1[0]
  for i in range(int(region1[1]), int(region1[2]), resolution):
    counterJ = 0
    for j in range(int(region2[1]), int(region2[2]), resolution):
      if((j-5*resolution < i and i < j+5*resolution) or (i-5*resolution < j and j < i+5*resolution)):
        pass
      else:
        try:
          value = float(matrixDict[":".join([chrom, str(i), str(j)])])
          matrix[counterI][counterJ] = matrix[counterI][counterJ] + value
        except Exception:
          try:
            value = float(matrixDict[":".join([chrom, str(j), str(i)])])
            matrix[counterI][counterJ] = matrix[counterI][counterJ] + value
          except Exception: pass
      counterJ += 1
    counterI += 1
  return 0

###################################################################################################
# Execution
###################################################################################################

# Reading all regions
regions = readRegions(regionFileName)

# Reading matrix dictionary
matrixDict = readMatrix(matrixFileName)

# Iterating on regions
counter = 0.0
matrix = [[0.0] * ((2*loopBins)+1) for e in range(0,((2*loopBins)+1))]
for r in regions:

  # Regions
  region1 = r[:3]
  region2 = r[3:]

  # Updating matrix
  updateMatrix(matrix, region1, region2, matrixDict, resolution)
 
  # Update average counter
  counter += 1.0

# Printing loop'O'gram
outputMatrixFile = open(outputMatrixFileName,"w")
for vec in matrix:
  toPrint = []
  for v in vec: toPrint.append(str(v / counter))
  outputMatrixFile.write("\t".join(toPrint)+"\n")
outputMatrixFile.close()


