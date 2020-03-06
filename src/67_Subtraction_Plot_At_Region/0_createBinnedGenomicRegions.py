
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
binSize = 5000000
chromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
outputFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/0_Regions/5MB.txt"

# Initialization
outputLocation = "/".join(outputFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def get_chrom_dict(chromSizesFileName):

  # Chromosome Dictionary
  chromSizesDict = dict()
  chromSizesFile = open(chromSizesFileName,"rU")
  for line in chromSizesFile:
    ll = line.strip().split("\t")
    chromSizesDict[ll[0]] = int(ll[1])
  chromSizesFile.close()
  chromList = sorted(chromSizesDict.keys())
  return chromList, chromSizesDict

def create_regions(binSize, chromSizesFileName, outputFileName):

  # Reading chrom dict
  chromList, chromSizesDict = get_chrom_dict(chromSizesFileName)

  # Main loop
  outputFile = open(outputFileName, "w")
  for chrom in chromList:
    for i in range(0, chromSizesDict[chrom]+binSize, binSize):
      p1 = i
      p2 = min(i + binSize, chromSizesDict[chrom])
      if(((p2 - p1) < binSize) or (p1 > chromSizesDict[chrom]) or (p2 > chromSizesDict[chrom])): continue
      outputFile.write("\t".join([chrom, str(p1), str(p2)])+"\n")
  outputFile.close()

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
create_regions(binSize, chromSizesFileName, outputFileName)


