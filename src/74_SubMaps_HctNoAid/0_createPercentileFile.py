
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from numpy import percentile

# Input
inputMatrixFileName = sys.argv[1]
outputPercFileName = sys.argv[2]

# Initialization
outputLocation = "/".join(outputPercFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Execution
###################################################################################################

# Read matrix into a list
matrixList = []
inputMatrixFile = open(inputMatrixFileName, "rU")
for line in inputMatrixFile:
  ll = line.strip().split("\t")
  matrixList.append(float(ll[3]))
inputMatrixFile.close()

# Calculate and write percentiles
outputPercFile = open(outputPercFileName, "w")
for i in range(0,101):
  perc = percentile(matrixList, i)
  outputPercFile.write("\t".join([str(e) for e in [i, perc]])+"\n")
outputPercFile.close()


