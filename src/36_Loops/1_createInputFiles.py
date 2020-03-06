
# Import
import os
import sys

# Input
resolution = int(sys.argv[1])
chromSizesFileName = sys.argv[2]
inputMatrixFileName = sys.argv[3]
outputFragmentFileName = sys.argv[4]
outputContactFileName = sys.argv[5]

# Initialization
chromosome = inputMatrixFileName.split("/")[-1].split(".")[0]

# Reading chrom sizes
chromSizesDict = dict()
chromSizesFile = open(chromSizesFileName,"rU")
for line in chromSizesFile:
  ll = line.strip().split("\t")
  chromSizesDict[ll[0]] = int(ll[1])
chromSizesFile.close()

# Reading input matrix
matrixDict = dict()
posDict = dict()
inputMatrixFile = open(inputMatrixFileName, "rU")
for line in inputMatrixFile:
  ll = line.strip().split("\t")
  key = ":".join([ll[0], ll[1], ll[2]])
  value = str(int(float(ll[3])))
  if(value == "0"): continue
  matrixDict[key] = value
  try: posDict[ll[1]] += 1
  except Exception: posDict[ll[1]] = 1
  try: posDict[ll[2]] += 1
  except Exception: posDict[ll[2]] = 1
inputMatrixFile.close()

# Create fragment file
outputFragmentFile = open(outputFragmentFileName, "w")
for i in range(resolution, chromSizesDict[chromosome]+resolution, resolution):
  if(i > chromSizesDict[chromosome]): continue
  mcc = "0"
  mapp = "0"
  try:
    mcc = str(posDict[str(i)])
    mapp = "1"
  except Exception: pass
  outputFragmentFile.write("\t".join([chromosome, "0", str(i), mcc, mapp])+"\n")
outputFragmentFile.close()

# Create matrix file
outputContactFile = open(outputContactFileName, "w")
for k in matrixDict.keys():
  kk = k.split(":")
  outputContactFile.write("\t".join([chromosome, kk[1], chromosome, kk[2], matrixDict[k]])+"\n")
outputContactFile.close()

# Compressing files
command = "gzip -f "+outputFragmentFileName
os.system(command)
command = "gzip -f "+outputContactFileName
os.system(command)


