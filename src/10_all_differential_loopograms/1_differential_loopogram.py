
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

# readRegionsSpecial
def readRegions(inputFileName):
  regions = []
  inputFile = open(inputFileName,"rU")
  for line in inputFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = ll[1]; p2 = ll[2]
    mid = (int(p1)+int(p2)) / 2.0
    r1 = str(int(math.floor(mid / float(resolution)) * float(resolution)) - ((loopBins/2) * int(resolution)))
    r2 = str(int(math.ceil(mid / float(resolution)) * float(resolution)) + ((loopBins/2) * int(resolution)))
    chrom = ll[3]; p1 = ll[4]; p2 = ll[5]
    mid = (int(p1)+int(p2)) / 2.0
    r3 = str(int(math.floor(mid / float(resolution)) * float(resolution)) - ((loopBins/2) * int(resolution)))
    r4 = str(int(math.ceil(mid / float(resolution)) * float(resolution)) + ((loopBins/2) * int(resolution)))
    region = [chrom, r1, r2, chrom, r3, r4]
    regions.append(region)
  inputFile.close()
  return regions

def updateMatrix(matrix, region1, region2, matrixDict, resolution):
  counterI = 0
  chrom = region1[0]
  for i in range(int(region1[1]),int(region1[2]),resolution):
    counterJ = 0
    for j in range(int(region2[1]),int(region2[2]),resolution):
      try:
        value = float(matrixDict[":".join([chrom,str(i),str(j)])])
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

# Reading matrix
matrixDict = dict()
matrixFile = open(matrixFileName,"rU")
for line in matrixFile:
  ll = line.strip().split("\t")
  matrixDict[":".join(ll[:3])] = float(ll[3])
matrixFile.close()

# Iterating on regions
counter = 0.0
matrix = [[0.0] * (loopBins+1) for e in range(0,(loopBins+1))]
for r in regions:

  # Regions
  region1 = r[:3]
  region2 = r[3:]

  # Updating matrix
  updateMatrix(matrix, region1, region2, matrixDict, resolution)
 
  # Update average counter
  counter += 1.0

  # Track
  print counter
  sys.stdout.flush()

# Printing loop'O'gram
outputMatrixFile = open(outputMatrixFileName,"w")
for vec in matrix:
  toPrint = []
  for v in vec: toPrint.append(str(v / counter))
  outputMatrixFile.write("\t".join(toPrint)+"\n")
outputMatrixFile.close()


