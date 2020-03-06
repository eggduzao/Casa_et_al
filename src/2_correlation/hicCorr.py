
# Import
import os
import sys
import math

# Input
expName = sys.argv[1]
juicerCommand = sys.argv[2]
kindOfMatrix = sys.argv[3]
kindOfNormalization = sys.argv[4]
unitOfResolution = sys.argv[5]
resolution = sys.argv[6]
smoothing = sys.argv[7]
maxDistInteract = sys.argv[8]
regionsFileName = sys.argv[9]
hicFileName1 = sys.argv[10]
hicFileName2 = sys.argv[11]
outputFileName = sys.argv[12]

# Functions

def extractMatrix(juicerCommand, kindOfMatrix, kindOfNormalization, hicFileName, region,
                  unitOfResolution, resolution, outFileName):
  command = " ".join([juicerCommand, "dump", kindOfMatrix, kindOfNormalization, hicFileName,
                      region, region, unitOfResolution, resolution, outFileName])
  os.system(command)

def extractMatrixException():
  

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

# Parameters
tempFileName = "./TEMP/"

# Iterating on regionsFileName
regionsFile = open(regionsFileName,"rU")
corrDict = dict()
for line in regionsFile:

  # Initialization
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]
  region = ":".join([chrom.split("chr")[-1],p1,p2])
  region2 = chrom+":"+p1+"-"+p2
  os.system("mkdir -p "+tempFileName)

  # Extracting matrix1
  sparseFileName1 = tempFileName+expName+"_"+region+"_sparse_1.txt"
  denseFileName1 = tempFileName+expName+"_"+region+"_dense_1.txt"
  try:
    extractMatrix(juicerCommand, kindOfMatrix, kindOfNormalization, hicFileName1, region,
                  unitOfResolution, resolution, sparseFileName1)
    sparseToDense(sparseFileName1, denseFileName1, chrom, p1, p2, resolution)
  except Exception: corrDict[region2] = "NA"

  # Extracting matrix2
  sparseFileName2 = tempFileName+expName+"_"+region+"_sparse_2.txt"
  denseFileName2 = tempFileName+expName+"_"+region+"_dense_2.txt"
  try:
    extractMatrix(juicerCommand, kindOfMatrix, kindOfNormalization, hicFileName2, region,
                  unitOfResolution, resolution, sparseFileName2)
    sparseToDense(sparseFileName2, denseFileName2, chrom, p1, p2, resolution)
  except Exception: corrDict[region2] = "NA"

  # Perform hicrep
  try:
    corrFileName = tempFileName+expName+"_"+region+"corr.txt"
    corrFileNameROut = tempFileName+expName+"_"+region+"corr.Rout"
    command = "R CMD BATCH '--args '"+denseFileName1+"' '"+denseFileName2+"' '"+corrFileName+"' '"+resolution+"' '"+smoothing+"' '"+maxDistInteract+" ./hicCorr.R "+corrFileNameROut
    os.system(command)
  except Exception: corrDict[region2] = "NA"

  # Read correlation
  try:
    corrFile = open(corrFileName,"rU")
    corrDict[region2] = corrFile.readline().strip()
    corrFile.close()
  except Exception: corrDict[region2] = "NA"

  # Termination
  os.system("rm -rf "+tempFileName)

regionsFile.close()

# Print correlations
outputFile = open(outputFileName,"w")
corrKeys = sorted(corrDict.keys())
for k in corrKeys:
  outputFile.write("\t".join([k,corrDict[k]])+"\n")
outputFile.close()


