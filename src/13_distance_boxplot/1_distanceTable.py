
# Import
import os
import sys

###################################################################################################
# Input
###################################################################################################

# Input
loopRegionsFileName = sys.argv[1]
outputFileName = sys.argv[2]

###################################################################################################
# Functions
###################################################################################################

def distance(r1, r2):
  x, y = sorted((r1, r2))
  if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
    return y[0] - x[1]
  return 0

###################################################################################################
# Execution
###################################################################################################

# Chrom list
loopRegionsFile = open(loopRegionsFileName, "rU")
outputFile = open(outputFileName, "w")
for line in loopRegionsFile:
  ll = line.strip().split("\t")
  r1 = [int(ll[1]), int(ll[2])]
  r2 = [int(ll[4]), int(ll[5])]
  outputFile.write(str(distance(r1, r2))+"\n")
loopRegionsFile.close()
outputFile.close()


