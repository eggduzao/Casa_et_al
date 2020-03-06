
# Import
import os
import sys

# Input
minV = int(sys.argv[1])
maxV = int(sys.argv[2])
inputFileName = sys.argv[3]
outputFileName = sys.argv[4]

# Execution
inputFile = open(inputFileName, "rU")
outputFile = open(outputFileName, "w")
inputFile.readline()
outputFile.write("\t".join(["START", "END", "ORDER"])+"\n")
for line in inputFile:
  ll = line.strip().split("\t")
  if(int(ll[0]) < minV or int(ll[1]) > maxV): continue
  outputFile.write("\t".join(ll)+"\n")


