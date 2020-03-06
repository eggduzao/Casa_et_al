
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
medThresh = sys.argv[3]
highThresh = sys.argv[4]
binsMiddle = int(sys.argv[5])
acceptedErrors = int(sys.argv[6])
minimumDistance = int(sys.argv[7])
thresholdFileName = sys.argv[8]
regionFileName = sys.argv[9]
matrixFileName = sys.argv[10]
outputRegionsFileName = sys.argv[11]
outputMatrixFileName = sys.argv[12]

###################################################################################################
# Functions
###################################################################################################

# getOverlap
# 0 = Intervals Overlap / Negative = int1 < int2 / Positive = int1 > int2
def getOverlap(a, b):
  if(min(a[1], b[1]) - max(a[0], b[0]) > 0): return 0
  else: return a[1] - b[0]

# readRegionsSpecial
def readRegions(inputFileName, chromList):
  regionDict = dict()
  originalRegionDict = dict()
  inputFile = open(inputFileName,"rU")
  for line in inputFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = ll[1]; p2 = ll[2]
    if(chrom not in chromList): continue
    mid = (int(p1)+int(p2)) / 2.0
    r1 = str(int(math.floor(mid / float(resolution)) * float(resolution)) - ((loopBins/2) * int(resolution)))
    r2 = str(int(math.ceil(mid / float(resolution)) * float(resolution)) + ((loopBins/2) * int(resolution)))
    region = [chrom, r1, r2]
    originalRegion = [chrom, p1, p2]
    try: regionDict[chrom].append(region)
    except Exception: regionDict[chrom] = [region]
    try: originalRegionDict[chrom].append(originalRegion)
    except Exception: originalRegionDict[chrom] = [originalRegion]
  inputFile.close()
  return regionDict, originalRegionDict

def updateMatrix(matrix, region1, region2, matrixDict, resolution, threshDict):
  counterI = 0
  chrom = region1[0]
  for i in range(int(region1[1]),int(region1[2]),resolution):
    counterJ = 0
    for j in range(int(region2[1]),int(region2[2]),resolution):
      try:
        value = float(matrixDict[":".join([chrom,str(i),str(j)])])
        if(value > threshDict["99"]): value = threshDict["99"]
        matrix[counterI][counterJ] = matrix[counterI][counterJ] + value
      except Exception: pass
      counterJ += 1
    counterI += 1
  return 0

def getMatrix(region1, region2, matrixDict, resolution, loopBins):
  counterI = 0
  chrom = region1[0]
  matrix = [[0.0] * (loopBins+1) for e in range(0,(loopBins+1))]
  for i in range(int(region1[1]),int(region1[2]),resolution):
    counterJ = 0
    for j in range(int(region2[1]),int(region2[2]),resolution):
      try: matrix[counterI][counterJ] = float(matrixDict[":".join([chrom,str(i),str(j)])])
      except Exception: matrix[counterI][counterJ] = 0.0
      counterJ += 1
    counterI += 1
  return matrix

###################################################################################################
# Execution
###################################################################################################

# Initializations
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Reading threshold
threshDict = dict()
thresholdFile = open(thresholdFileName, "rU")
for line in thresholdFile:
  ll = line.strip().split("\t")
  threshDict[ll[0]] = float(ll[1])
thresholdFile.close()

# Reading matrix
matrixDict = dict()
matrixFile = open(matrixFileName,"rU")
for line in matrixFile:
  ll = line.strip().split("\t")
  matrixDict[":".join(ll[:3])] = float(ll[3])
matrixFile.close()

# Reading all regions
regionDict, originalRegionDict = readRegions(regionFileName, chrList)

# Iterating on regions
counter = 0.0
matrix = [[0.0] * (loopBins+1) for e in range(0,(loopBins+1))]
outputRegionsFile = open(outputRegionsFileName, "w")
for chrom in sorted(regionDict.keys()):
  
  # Reading all regions
  regionList = regionDict[chrom]
  originalRegionList = originalRegionDict[chrom]

  # Iterating on regions
  for i in range(0,len(regionList)-1):
    for j in range(i+1,len(regionList)):

      # Regions
      region1 = regionList[i]
      region2 = regionList[j]

      # Minimum distance check
      int1 = [int(e) for e in region1[1:]]; int2 = [int(e) for e in region2[1:]]
      ovl = abs(getOverlap(int1, int2))
      if(ovl < minimumDistance): continue

      # Signal check
      midValue = 0
      m1 = ((int1[0]+int1[1])/2)
      m2 = ((int2[0]+int2[1])/2)
      hr = resolution/2
      r = resolution
      currErrors = 0
      for ii in range(m1-hr-(binsMiddle*r), m1-hr+(binsMiddle*r)+1, r):
        for jj in range(m2-hr-(binsMiddle*r), m2-hr+(binsMiddle*r)+1, r):
          key = ":".join([chrom,str(ii),str(jj)])
          try:
            if((ii == (m1-hr)) and (jj == (m2-hr))):
              midValue = matrixDict[key]
              if(midValue < threshDict[highThresh]): currErrors = acceptedErrors + 1
            else:
              if(midValue < threshDict[medThresh]): currErrors += 1
          except Exception: currErrors += 1
          if(currErrors > acceptedErrors): break
        if(currErrors > acceptedErrors): break
      if(currErrors > acceptedErrors): continue

      # Updating matrix
      updateMatrix(matrix, region1, region2, matrixDict, resolution, threshDict)
 
      # Update average counter
      counter += 1.0

      # Writing region
      originalRegion1 = [str(e) for e in originalRegionList[i]]
      originalRegion2 = [str(e) for e in originalRegionList[j]]
      outputRegionsFile.write("\t".join(originalRegion1+originalRegion2+[str(midValue)])+"\n")

      # Track
      print chrom, i, j
      sys.stdout.flush()

outputRegionsFile.close()

# Printing loop'O'gram
outputMatrixFile = open(outputMatrixFileName,"w")
for vec in matrix:
  toPrint = []
  for v in vec: toPrint.append(str(v / counter))
  outputMatrixFile.write("\t".join(toPrint)+"\n")
outputMatrixFile.close()


