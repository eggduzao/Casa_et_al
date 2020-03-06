
# Import
import os
import sys
from glob import glob

# Fixed parameters
juicerCommand = "juicertools"
unitOfResolution = "BP"
regionsFileName = "/home/egg/Projects/hic_corr/data/testRegions.txt"
hicFileName1 = "/home/egg/Projects/hic_corr/data/yulia_data/79643_30.hic"
hicFileName2 = "/home/egg/Projects/hic_corr/data/yulia_data/79644_30.hic"
testLocation = "/home/egg/Projects/hic_corr/TEST/"
outFileName = testLocation+"results"

# Variable parameters
kindOfMatrix = ["observed", "oe"]
kindOfNormalization = ["NONE", "VC", "VC_SQRT", "KR"]
resolution = ["250000", "100000", "50000", "25000", "10000", "5000"]
smoothing = ["0", "1", "2"]
maxDistInteract = ["1000000", "2500000", "5000000"]

# Execution
for a in kindOfMatrix:
  for b in kindOfNormalization:
    for c in resolution:
      for d in smoothing:
        for e in maxDistInteract:
          expName = "_".join([a,b,c,d,e])
          outputFileName = testLocation+expName+".info"
          command = "./run.sh "+" ".join([expName, juicerCommand, a, b, unitOfResolution,
                    c, d, e, regionsFileName, hicFileName1, hicFileName2, outputFileName])
          os.system(command)

# Creating final list
inList = sorted(glob(testLocation+"*.info"))
outFileUnsortName = outFileName+"unsorted.info"
outFileUnsort = open(outFileUnsortName,"w")
for inFileName in inList:
  name = inFileName.split("/")[-1].split(".")[0]
  inFile = open(inFileName,"rU")
  for line in inFile:
    ll = line.strip().split("\t")
    value = ll[1]+"\t"+ll[2]
    break
  outFileUnsort.write("\t".join([name,value])+"\n")
  inFile.close()
outFileUnsort.close()
  
# Sorting final list
outFileSortName = outFileName+".txt"
os.system("sort -k2,2rn "+outFileUnsortName+" > "+outFileSortName)

# Termination
# os.system("rm *.info")


