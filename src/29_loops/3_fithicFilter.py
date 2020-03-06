
# Import
import os
import sys

# Input
qValueThreshold = float(sys.argv[1])
inputContactFileName = sys.argv[2]
tempLoc = sys.argv[3]
outputContactFileName = sys.argv[4]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

# Creating temporary uncompressed file
uncompressedFileName = tempLoc + "uncompressedFileName.txt"
command = "gzip -cd "+inputContactFileName+" > "+uncompressedFileName
os.system(command)

# Thresholding file
inputContactFile = open(uncompressedFileName, "rU")
tempContactFileName = tempLoc + "tempContactFileName.txt"
tempContactFile = open(tempContactFileName, "w")
inputContactFile.readline()
for line in inputContactFile:
  ll = line.strip().split("\t")
  qvalue = float(ll[6])
  if(qvalue > qValueThreshold): continue
  tempContactFile.write("\t".join([ll[0], ll[1], ll[3], ll[4]])+"\n")
inputContactFile.close()
tempContactFile.close()

# Sorting output
command = "sort -k1,1 -k2,2n "+tempContactFileName+" > "+outputContactFileName
os.system(command)

# Removing temp
command = "rm -rf "+tempLoc
os.system(command)


