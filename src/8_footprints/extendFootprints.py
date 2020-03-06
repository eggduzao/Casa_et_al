
#Import
import os
import sys

# Input
ext = "25"
inputFpList = ["xxx"]
outputFpList = ["xxx"]

# Initialization
for i in range(0,len(inputFpList)):
  inputFpFileName = inputFpList[i]
  outputFpFileName = outputFpList[i]
  inputFpFile = open(inputFpFileName, "rU")
  outputFpFile = open(outputFpFileName, "rU")
  for line in inputFpFile:
    ll = line.strip().split("\t")
    mid = (int(ll[1])+int(ll[2]))/2
    outputFpFile.write("\t".join([ll[0], mid-ext, mid+ext])+"\n")
  inputFpFile.close()
  outputFpFile.close()


