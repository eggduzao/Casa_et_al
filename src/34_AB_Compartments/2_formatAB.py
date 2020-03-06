
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

# Input
chromosome = sys.argv[1]
resolution = int(sys.argv[2])
chromSizesFileName = sys.argv[3]
enhChromSizesFileName = sys.argv[4]
eigenFileName = sys.argv[5]
dnaseFileName = sys.argv[6]
tempLocation = sys.argv[7]
outBedFileName = sys.argv[8]
outBwFileName = sys.argv[9]
outCountFileName = sys.argv[10]

# Initialization
pseudocount = 0.0001
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def fetchTotalReadsBam(bamFile, region):
  returnN = 0
  try:
    for read in bamFile.fetch(region[0], region[1], region[2]): returnN += 1
  except Exception: pass
  return returnN

###################################################################################################
# Basic Structures / Correlation with DNase to check signal
###################################################################################################

# Chromosome Dictionary
chromSizesDict = dict() # CHROMOSOME -> SIZE
chromSizesFile = open(chromSizesFileName,"rU")
for line in chromSizesFile:
  ll = line.strip().split("\t")
  chromSizesDict[ll[0]] = int(ll[1])
chromSizesFile.close()

# Creating eigenvector list
eigenList = [] # [CHR:p1:p2, value]
posReadCount = 0
posCount = 0
negReadCount = 0
negCount = 0
position = 0
minNeg = 999
maxPos = -999
eigenFile = open(eigenFileName,"rU")
dnaseFile = Samfile(dnaseFileName, "rb")
for line in eigenFile:
  ll = line.strip().split("\t")
  p1 = position; p2 = min(position+resolution, chromSizesDict[chromosome])
  key = ":".join([chromosome, str(p1), str(p2)])
  if(ll[0] == "NaN" or ll[0] == "Inf" or ll[0] == "-Inf"): ll[0] = 0
  value = float(ll[0])
  if(value < 0): value -= pseudocount
  if(value > 0): value += pseudocount
  eigenVec = [key, value]
  if(value > 0):
    posReadCount += float(fetchTotalReadsBam(dnaseFile, [chromosome, p1, p2]))
    posCount += 1
    if(value > maxPos): maxPos = value
  elif(value < 0):
    negReadCount += float(fetchTotalReadsBam(dnaseFile, [chromosome, p1, p2]))
    negCount += 1
    if(value < minNeg): minNeg = value
  position += resolution
  eigenList.append(eigenVec)
eigenFile.close()
dnaseFile.close()
posReadCount = posReadCount/posCount
negReadCount = negReadCount/negCount
maxPos = round(maxPos * 100,4)
minNeg = round(minNeg * 100,4)

# Checking if signal change is needed
if(negReadCount > posReadCount):
  for i in range(0,len(eigenList)):
    if(eigenList[i][1] == 0):
      eigenList[i][1] = -maxPos
      continue
    eigenList[i][1] = round(-eigenList[i][1] * 100,4)
else:
  for i in range(0,len(eigenList)):
    if(eigenList[i][1] == 0):
      eigenList[i][1] = minNeg
      continue
    eigenList[i][1] = round(eigenList[i][1] * 100,4)

###################################################################################################
# Writing output files
###################################################################################################

# Writing counts
outCountFile = open(outCountFileName, "w")
outCountFile.write("\t".join(["POSITIVE_COUNTS", "NEGATIVE_COUNTS"])+"\n")
outCountFile.write("\t".join([str(round(e,2)) for e in [posReadCount, negReadCount]])+"\n")
outCountFile.close()

# Writing A and B compartments
outBedFile = open(outBedFileName, "w")
for i in range(0,len(eigenList)):
  kk = eigenList[i][0].split(":")
  chrom = kk[0]; p1 = kk[1]; p2 = kk[2]
  name = "NA"
  if(eigenList[i][1] > 0): name = "A"
  elif(eigenList[i][1] < 0): name = "B"
  value = str(eigenList[i][1])
  outBedFile.write("\t".join([chrom, p1, p2, name, value, "."])+"\n")
outBedFile.close()

# Writing bw file
tempWigFileName = tempLocation+"tempWigFile.wig"
tempWigFile = open(tempWigFileName, "w")
tempWigFile.write("fixedStep chrom="+chromosome+" start=1 step="+str(resolution)+" span="+str(resolution-1)+"\n")
for i in range(0,len(eigenList)):
  tempWigFile.write(str(eigenList[i][1])+"\n")
tempWigFile.close()
command = "wigToBigWig "+" ".join([tempWigFileName, enhChromSizesFileName, outBwFileName])
os.system(command)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


