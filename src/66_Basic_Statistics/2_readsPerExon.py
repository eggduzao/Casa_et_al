
# Import
import os
import sys
from pysam import Samfile

###################################################################################################
# INPUT
###################################################################################################

# Input
regionFileName = sys.argv[1]
bamFileName = sys.argv[2]
outputFileName = sys.argv[3]

# Create output location
outputLocation = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# FUNCTIONS
###################################################################################################

def fetchTotalReadsBam(bamFile, region):
  returnN = 0
  for read in bamFile.fetch(region[0], region[1], region[2]): returnN += 1
  return returnN

###################################################################################################
# EXECUTION
###################################################################################################

summ = 0.0; total = 0.0
bamFile = Samfile(bamFileName, "rb")
regionFile = open(regionFileName, "rU")
outputFile = open(outputFileName, "w")
for line in regionFile:
  ll = line.strip().split("\t")
  lll = ll[3].split(":")
  if(lll[0] != "EXON" or lll[2] != "ACTIVE"): continue
  try:
    summ += fetchTotalReadsBam(bamFile, [ll[0], int(ll[1]), int(ll[2])])
  except Exception: continue
  total += 1.0
bamFile.close()
regionFile.close()
outputFile.write("\t".join(["SUMM", "TOTAL", "FRACTION"])+"\n")
outputFile.write("\t".join([str(summ), str(total), str(round(summ/total,2))])+"\n")
outputFile.close()


