
# Import
import os
import sys

# Input


# Functions

def extractMatrix(juicerCommand, kindOfMatrix, kindOfNormalization, hicFileName, region,
                  unitOfResolution, resolution, outFileName):
  command = " ".join([juicerCommand, "dump", kindOfMatrix, kindOfNormalization, hicFileName,
                      region, region, unitOfResolution, resolution, outFileName])
  os.system(command)

def sparseToDense(sparseFileName, denseFileName, chrom, p1, p2, resolution):
  helpDict = dict()
  k1 = int(p1); k2 = int(p2); res = int(resolution)
  sparseFile = open(sparseFileName,"rU")
  for line in sparseFile:
    ll = line.strip().split("\t")
    helpDict[":".join([ll[0],ll[1]])] = ll[2]
  sparseFile.close()
  denseFile = open(denseFileName,"w")
  for i in range(k1, k2, res):
    toWrite = [chrom, i, i+res]
    for j in range(k1, k2, res):
      minP = str(min(i, j))
      maxP = str(max(i, j))
      try: value = float(helpDict[":".join([minP,maxP])])
      except Exception: value = 0.0
      if(math.isnan(value) or math.isinf(value)): value = 0
      toWrite.append(value)
    denseFile.write("\t".join([str(e) for e in toWrite])+"\n")   
  denseFile.close()

def recreteHic(): pass

# Input
expName = sys.argv[1]
juicerCommand = sys.argv[2]
kindOfMatrix = sys.argv[3]
kindOfNormalization = sys.argv[4]
unitOfResolution = sys.argv[5]
resolution = sys.argv[6]
regionsFileName = sys.argv[7]
hicFileName1 = sys.argv[8]
hicFileName2 = sys.argv[9]
outputFileName = sys.argv[10]

# Parameters
os.system("mkdir -p "+tempFileName)

# Iterating on regionsFileName
regionsFile = open(regionsFileName,"rU")
for line in regionsFile:

  # Initialization
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]
  region = ":".join([chrom.split("chr")[-1],p1,p2])

  # Extracting matrix1
  try:
    sparseFileName1 = tempFileName+expName+"_"+region+"_sparse_1.txt"
    denseFileName1 = tempFileName+expName+"_"+region+"_dense_1.txt"
    extractMatrix(juicerCommand, kindOfMatrix, kindOfNormalization, hicFileName1, region,
                  unitOfResolution, resolution, sparseFileName1)
    sparseToDense(sparseFileName1, denseFileName1, chrom, p1, p2, resolution)
  except Exception: pass

  # Extracting matrix2
  try:
    sparseFileName2 = tempFileName+expName+"_"+region+"_sparse_2.txt"
    denseFileName2 = tempFileName+expName+"_"+region+"_dense_2.txt"
    extractMatrix(juicerCommand, kindOfMatrix, kindOfNormalization, hicFileName2, region,
                  unitOfResolution, resolution, sparseFileName2)
    sparseToDense(sparseFileName2, denseFileName2, chrom, p1, p2, resolution)
  except Exception: pass

regionsFile.close()




# Termination
os.system("rm -rf "+tempFileName)


