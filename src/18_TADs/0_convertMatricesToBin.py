
# Import
import os
import sys
from pysam import Samfile

###################################################################################################
# INPUT
###################################################################################################

# Input
chromosome = sys.argv[1]
resolution = int(sys.argv[2])
matrixFileName = sys.argv[3]
chromSizesFileName = sys.argv[4]
outputFileName = sys.argv[5]

###################################################################################################
# EXECUTION
###################################################################################################

# Reading chrom sizes
chromSizesDict = dict()
chromSizesFile = open(chromSizesFileName,"rU")
for line in chromSizesFile:
  ll = line.strip().split("\t")
  chromSizesDict[ll[0]] = int(ll[1])
chromSizesFile.close()
chromList = chromSizesDict.keys()

# Creating bin dictionary
binDict = dict()
counter = 0
for i in range(0,chromSizesDict[chromosome],resolution):
  binDict[str(i)] = str(counter)
  counter += 1

# Writing output
matrixFile = open(matrixFileName, "rU")
outputFile = open(outputFileName, "w")
for line in matrixFile:
  ll = line.strip().split("\t")
  outputFile.write("\t".join([ll[0], binDict[ll[1]], binDict[ll[2]], ll[3]])+"\n")
  
# Closing files
matrixFile.close()
outputFile.close()


