
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
expThresh = float(sys.argv[1])
expressionFileName = sys.argv[2]
outputActiveFileName = sys.argv[3]
outputInactiveFileName = sys.argv[4]
  
###################################################################################################
# Execution
###################################################################################################

# Thresholding expression file
expressionFile = open(expressionFileName, "rU")
outputActiveFile = open(outputActiveFileName, "w")
outputInactiveFile = open(outputInactiveFileName, "w")
for line in expressionFile:
  ll = line.strip().split("\t")
  gene = ll[0]; exp = float(ll[1])
  if(exp <= expThresh): outputInactiveFile.write("\t".join(ll)+"\n")
  else: outputActiveFile.write("\t".join(ll)+"\n")
expressionFile.close()
outputActiveFile.close()
outputInactiveFile.close()


