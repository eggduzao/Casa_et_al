
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import math
import pyBigWig
import numpy as np
from pysam import Samfile

# Input
plotExt = int(sys.argv[1])
signalExtU = int(sys.argv[2])
signalExtD = int(sys.argv[3])
bedFileName = sys.argv[4]
bamFileName = sys.argv[5]
outputFileName = sys.argv[6]

# Initialization
outLoc = "/".join(outputFileName.split("/")[:-1])+"/"
command = "mkdir -p "+outLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def fetchSignalBam(bamFile, region, signalExtU, signalExtD):
  regionLen = (region[2] - region[1])
  returnVec = [0.0] * regionLen
  for read in bamFile.fetch(region[0], region[1], region[2]):
    readLoc = read.reference_start - region[1]
    if(not read.is_reverse):
      for i in range(max(readLoc-signalExtU, 0), min(readLoc+signalExtD, len(returnVec))): returnVec[i] += 1.0
    else:
      for i in range(max(readLoc-signalExtD, 0), min(readLoc+signalExtU, len(returnVec))): returnVec[i] += 1.0
  return returnVec

###################################################################################################
# Creating table
###################################################################################################

# Initialization
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]
signalFile = Samfile(bamFileName, "rb")

# Fetching the bam signal in all categories
# REGION_CHR, REGION_P1, REGION_P2, [SIGNAL...]
bedFile = open(bedFileName, "rU")
outputFile = open(outputFileName,"w")
for line in bedFile:

  # Initialization
  ll = line.strip().split("\t")
  if(ll[0] not in chrList): continue
  chrom1 = ll[0]; px1 = int(ll[1]); py1 = int(ll[2])
  chrom2 = ll[3]; px2 = int(ll[4]); py2 = int(ll[5])
  vector1 = [chrom1, px1, py1]
  vector2 = [chrom2, px2, py2]

  # Fetching signal 1
  mid = (px1+py1)/2
  region = [chrom1, mid-plotExt, mid+plotExt]
  signal = fetchSignalBam(signalFile, region, signalExtU, signalExtD)
  vector1 = vector1 + signal
  outputFile.write("\t".join([str(e) for e in vector1])+"\n")

  # Fetching signal 2
  mid = (px2+py2)/2
  region = [chrom2, mid-plotExt, mid+plotExt]
  signal = fetchSignalBam(signalFile, region, signalExtU, signalExtD)
  vector2 = vector2 + signal
  outputFile.write("\t".join([str(e) for e in vector2])+"\n")

# Closing all files
bedFile.close()
outputFile.close()
signalFile.close()


