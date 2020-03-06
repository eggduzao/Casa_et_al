###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
chromosome = sys.argv[1]
meanHalfSize = int(sys.argv[2])
resolution = int(sys.argv[3])
chromSizesFileName = sys.argv[4]
matrixFileName = sys.argv[5]
outputFileName = sys.argv[6]

# Creating output folder
outLoc = "/".join(outputFileName.split("/")[:-1])+"/"
os.system("mkdir -p "+outLoc)

###################################################################################################
# Execution
###################################################################################################

# Reading chrom sizes
chromSizesDict = dict()
chromSizesFile = open(chromSizesFileName,"rU")
for line in chromSizesFile:
  ll = line.strip().split("\t")
  chromSizesDict[ll[0]] = int(ll[1]) - (int(ll[1])%resolution)
chromSizesFile.close()
chromList = chromSizesDict.keys()

# Reading matrix
matrixDict = dict()
matrixFile = open(matrixFileName,"rU")
for line in matrixFile:
  ll = line.strip().split("\t")
  matrixDict[":".join(ll[:3])] = float(ll[3])
matrixFile.close()

# Iterating on all possible genomic positions
resDict = dict()
for i in range(0, chromSizesDict[chromosome]-resolution, resolution):
  for j in range(i, chromSizesDict[chromosome]+1, resolution):
    distance = j - i
    try: value = matrixDict[":".join([chromosome, str(i), str(j)])]
    except Exception: value = 0.0
    try: resDict[distance].append(value)
    except Exception: resDict[distance] = [value]
sortedDistances = sorted(resDict.keys())

# Updating resDict
resDictNorm = dict()
outFile = open(outputFileName,"w")
for k in sortedDistances:
  value = sum(resDict[k])
  value2 = sum(resDict[k]) / len(resDict[k])
  resDict[k] = value
  resDictNorm[k] = value
  outFile.write("\t".join([str(k), str(value), str(value2)])+"\n")
outFile.close()


