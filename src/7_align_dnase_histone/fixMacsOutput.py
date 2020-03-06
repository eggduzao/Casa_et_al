
#Import
import os
import sys
import numpy as np

###################################################################################################
# Input
###################################################################################################

# Input
peakFileName = sys.argv[1]
summitFileName = sys.argv[2]

# Parameters
percentileToFilter = 25

###################################################################################################
# Functions
###################################################################################################

def getPercentile(percentileToFilter, fileName, loc):
  inFile = open(fileName,"rU")
  vec = []
  for line in inFile:
    ll = line.strip().split("\t")
    vec.append(float(ll[loc]))
  inFile.close()
  vec = np.array(vec)
  return np.percentile(vec, percentileToFilter)

###################################################################################################
# Execution
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Fixing peaks
pName = peakFileName.split(".")[0]
outputPeakFileNameT = pName+"_filterT.narrowPeak"
outputPeakFileName = pName+"_filter.narrowPeak"
command = "mergeBed -c 4,5,6,7,8,9,10 -o collapse,mean,distinct,mean,mean,mean,median -i "+peakFileName+" > "+outputPeakFileNameT
os.system(command)
loc = 7
minScore = getPercentile(percentileToFilter, outputPeakFileNameT, loc)
inFile = open(outputPeakFileNameT, "rU")
outFile = open(outputPeakFileName, "w")
for line in inFile:
  ll = line.strip().split("\t")
  if(ll[0] not in chrList): continue
  if(ll[loc] < minScore): continue
  outFile.write(line)
inFile.close()
outFile.close()
os.system("rm "+outputPeakFileNameT)

# Fixing summits
sName = summitFileName.split(".")[0]
outputSummitFileName = sName+"_filter.bed"
command = "intersectBed -wa -u -a "+summitFileName+" -b "+outputPeakFileName+" > "+outputSummitFileName
os.system(command)


