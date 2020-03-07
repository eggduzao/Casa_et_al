
# Import
import os
import sys

###################################################################################################
# INPUT
###################################################################################################

# Input
chromosome = sys.argv[1]
inputFileName = sys.argv[2]
temporaryLocation = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
command = "mkdir -p "+temporaryLocation
os.system(command)
outputLocation = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# EXECUTION
###################################################################################################

# Fetching chromosome only
inputFile = open(inputFileName, "rU")
tempUnsortedFileName = temporaryLocation + "tempUnsortedFileName.txt"
tempUnsortedFile = open(tempUnsortedFileName, "w")
for line in inputFile:
  ll = line.strip().split("\t")
  if(ll[0] != chromosome): continue
  tempUnsortedFile.write(line)
inputFile.close()
tempUnsortedFile.close()

# Sorting file
command = "sort -k1,1 -k2,2n -k3,3n " + tempUnsortedFileName + " > " + outputFileName
os.system(command)


