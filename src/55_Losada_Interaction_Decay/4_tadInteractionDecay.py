###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
chromosome = sys.argv[1]
resolution = int(sys.argv[2])
chromSizesFileName = sys.argv[3]
matrixMinusFileName = sys.argv[4]
matrixPlusFileName = sys.argv[5]
outputFileName = sys.argv[6]

# Creating output folder
maxDistance = 2000000
outLoc = "/".join(outputFileName.split("/")[:-1])+"/"
os.system("mkdir -p "+outLoc)

###################################################################################################
# Functions
###################################################################################################

def standardizeDict(d):
  minV = 99999999.
  maxV = -99999999.
  for k in d.keys():
    if(d[k] > maxV): maxV = d[k]
    if(d[k] < minV): minV = d[k]
  newDict = dict()
  for k in d.keys(): newDict[k] = (d[k] - minV) / (maxV - minV)
  return newDict

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

# Reading minus matrix
matrixMinusDict = dict()
matrixMinusFile = open(matrixMinusFileName,"rU")
for line in matrixMinusFile:
  ll = line.strip().split("\t")
  matrixMinusDict[":".join(ll[:3])] = float(ll[3])
matrixMinusFile.close()

matrixMinusDict = standardizeDict(matrixMinusDict)

# Reading plus matrix
matrixPlusDict = dict()
matrixPlusFile = open(matrixPlusFileName,"rU")
for line in matrixPlusFile:
  ll = line.strip().split("\t")
  matrixPlusDict[":".join(ll[:3])] = float(ll[3])
matrixPlusFile.close()

matrixPlusDict = standardizeDict(matrixPlusDict)

# Creating fold change interaction-decay
fcDict = dict()
allKeys = list(set(matrixMinusDict.keys()).union(set(matrixPlusDict.keys())))
for k in allKeys:
  try: mValue = matrixMinusDict[k]
  except Exception:
    fcDict[k] = 0.0
    continue
  try: pValue = matrixPlusDict[k]
  except Exception:
    fcDict[k] = 0.0
    continue
  if(mValue > 0): fcDict[k] = pValue / mValue
  else: continue

# Iterating on all possible genomic positions (until 2 Mb)
resDict = dict()
for i in range(0, chromSizesDict[chromosome]-maxDistance, resolution):
  for j in range(i, i+maxDistance, resolution):
    distance = (j - i)/resolution
    try: value = fcDict[":".join([chromosome, str(i), str(j)])]
    except Exception: value = 0.0
    try: resDict[distance].append(value)
    except Exception: resDict[distance] = [value]
sortedDistances = sorted(resDict.keys())

# Updating resDict
outFile = open(outputFileName,"w")
outFile.write("\t".join(["BIN", "TOTAL_VALUE", "NORM_VALUE"])+"\n")
for k in sortedDistances:
  value = sum(resDict[k])
  value2 = sum(resDict[k]) / len(resDict[k])
  outFile.write("\t".join([str(k), str(value), str(value2)])+"\n")
outFile.close()


