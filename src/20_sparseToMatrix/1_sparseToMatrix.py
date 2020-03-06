
# Import
import os
import sys

###################################################################################################
# INPUT
###################################################################################################

# Input
start = int(sys.argv[1])
end = int(sys.argv[2])
resolution = int(sys.argv[3])
sparseFileName = sys.argv[4]
outputFileName = sys.argv[5]

###################################################################################################
# EXECUTION
###################################################################################################

# Creating output matrix
nBins = (end - start) / resolution
matrix = [["0.0"] * nBins for e in range(0,nBins)]

# Populating matrix
sparseFile = open(sparseFileName, "rU")
for line in sparseFile:
  ll = line.strip().split("\t")
  b1 = (int(ll[1])-start)/resolution
  b2 = (int(ll[2])-start)/resolution
  if(b1 < 0 or b1 >= len(matrix) or b2 < 0 or b2 >= len(matrix)): continue
  matrix[b1][b2] = str(round(float(ll[3]),4))
  matrix[b2][b1] = str(round(float(ll[3]),4))

# Writing matrix
outputFile = open(outputFileName, "w")
for i in range(0,nBins):
  vec = []
  for j in range(0,nBins):
    vec.append(matrix[i][j])
  outputFile.write("\t".join(vec)+"\n")


