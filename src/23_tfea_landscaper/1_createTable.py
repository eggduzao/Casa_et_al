
# Import
import os
import sys

# Input
inputFileName1 = sys.argv[1]
inputFileName2 = sys.argv[2]
outputFileName = sys.argv[3]

# Reading input files
inDict1 = dict()
inputFile = open(inputFileName1, "rU")
inputFile.readline()
for line in inputFile:
  ll = line.strip().split(",")
  inDict1[ll[0]] = float(ll[2])
inputFile.close()
inDict2 = dict()
inputFile = open(inputFileName2, "rU")
inputFile.readline()
for line in inputFile:
  ll = line.strip().split(",")
  inDict2[ll[0]] = float(ll[2])
inputFile.close()

# Writing output
outputFile = open(outputFileName, "w")
for k in sorted(inDict1.keys()):
  outputFile.write("\t".join([k, str(inDict1[k] - inDict2[k])])+"\n")
outputFile.close()


