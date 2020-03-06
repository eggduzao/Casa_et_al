
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import numpy as np

# Input
plusAuxinFileName = sys.argv[1]
minusAuxinFileName = sys.argv[2]
specificFileName = sys.argv[3]
sensitiveFileName = sys.argv[4]
bothFileName = sys.argv[5]

###################################################################################################
# Execution
###################################################################################################

# Plus dict
plusDict = dict()
plusAuxinFile = open(plusAuxinFileName, "rU")
for line in plusAuxinFile:
  ll = line.strip().split("\t")
  plusDict[":".join(ll[:-1])] = ll[-1]
plusAuxinFile.close()

# Minus dict
minusDict = dict()
minusAuxinFile = open(minusAuxinFileName, "rU")
for line in minusAuxinFile:
  ll = line.strip().split("\t")
  minusDict[":".join(ll[:-1])] = ll[-1]
minusAuxinFile.close()

# Differential
plusKeys = sorted(plusDict.keys())
minusKeys = sorted(minusDict.keys())
specificKeys = list(np.setdiff1d(plusKeys, minusKeys, assume_unique=True))
sensitiveKeys = list(np.setdiff1d(minusKeys, plusKeys, assume_unique=True))
bothKeys = list(set(minusKeys).intersection(plusKeys))

# Writing to temp file
tempSpecificFileName = "./tempSpecificFileName.txt"
tempSensitiveFileName = "./tempSensitiveFileName.txt"
tempBothFileName = "./tempBothFileName.txt"
tempSpecificFile = open(tempSpecificFileName, "w")
tempSensitiveFile = open(tempSensitiveFileName, "w")
tempBothFile = open(tempBothFileName, "w")
for k in specificKeys: tempSpecificFile.write("\t".join(k.split(":")+[plusDict[k]])+"\n")
for k in sensitiveKeys: tempSensitiveFile.write("\t".join(k.split(":")+[minusDict[k]])+"\n")
for k in bothKeys: tempBothFile.write("\t".join(k.split(":")+[minusDict[k]])+"\n")
tempSpecificFile.close()
tempSensitiveFile.close()
tempBothFile.close()

# Sort based on score column
os.system("sort -k7,7gr "+tempSpecificFileName+" > "+specificFileName)
os.system("sort -k7,7gr "+tempSensitiveFileName+" > "+sensitiveFileName)
os.system("sort -k7,7gr "+tempBothFileName+" > "+bothFileName)

# Remove temporary files
os.system("rm "+tempSpecificFileName+" "+tempSensitiveFileName+" "+tempBothFileName)


