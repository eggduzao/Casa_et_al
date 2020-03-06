
#Import
import os
import sys

# Input
inputNameList = sys.argv[1].split(",")
inputFileNameList = sys.argv[2].split(",")
outputFileName = sys.argv[3]

# Create temp and output location
outputLocation = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

# Iterating on input files to fetch distribution
distMat = [[] for e in range(0, len(inputFileNameList))]
for i in range(0, len(inputFileNameList)):
  inputFile = open(inputFileNameList[i], "rU")
  for line in inputFile:
    ll = line.strip().split("\t")
    distMat[i].append(str(int(ll[2]) - int(ll[1])))
  inputFile.close()

# Writing output
maxL = max([len(e) for e in distMat])
outputFile = open(outputFileName, "w")
outputFile.write("\t".join(inputNameList)+"\n")
for j in range(0, maxL):
  vec = []
  for i in range(0, len(distMat)):
    try: vec.append(distMat[i][j])
    except Exception: vec.append("NA")
  outputFile.write("\t".join(vec)+"\n")
outputFile.close()


