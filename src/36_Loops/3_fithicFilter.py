
# Import
import os
import sys
from numpy import array, percentile

# Input
pValueThreshold = float(sys.argv[1])
percentileThreshold = float(sys.argv[2])
inputContactFileName = sys.argv[3]
tempLoc = sys.argv[4]
outputContactFileName = sys.argv[5]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

# Creating temporary uncompressed file
uncompressedFileName = tempLoc + "uncompressedFileName.txt"
command = "gzip -cd "+inputContactFileName+" > "+uncompressedFileName
os.system(command)

# Thresholding by p-value
inputContactFile = open(uncompressedFileName, "rU")
inputContactFile.readline()
signalVec = []
for line in inputContactFile:
  ll = line.strip().split("\t")
  qvalue = float(ll[6])
  if(qvalue > pValueThreshold): continue
  signalVec.append(float(ll[4]))
inputContactFile.close()

# Calculating percentile
signalVec = array(signalVec)
percentileValue = percentile(signalVec, percentileThreshold)
print percentileValue

# Thresholding file
inputContactFile = open(uncompressedFileName, "rU")
tempContactFileName = tempLoc + "tempContactFileName.txt"
tempContactFile = open(tempContactFileName, "w")
inputContactFile.readline()
for line in inputContactFile:
  ll = line.strip().split("\t")
  qvalue = float(ll[6])
  signal = float(ll[4])
  if((qvalue > pValueThreshold) or (signal < percentileValue)): continue
  tempContactFile.write("\t".join([ll[0], ll[1], ll[3], ll[4]])+"\n")
inputContactFile.close()
tempContactFile.close()

# Sorting output
command = "sort -k1,1 -k2,2n "+tempContactFileName+" > "+outputContactFileName
os.system(command)

# Removing temp
command = "rm -rf "+tempLoc
os.system(command)


