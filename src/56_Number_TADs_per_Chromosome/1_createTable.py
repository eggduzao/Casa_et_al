
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from glob import glob

# Input
prefix = sys.argv[1]
categoryList = sys.argv[2].split(",")
inputLocation = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
outputLocation = "/".join(outputFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def file_len(fname):
  i = -1
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i + 1

###################################################################################################
# Execution
###################################################################################################

# Creating table
table = [[] for e in categoryList]
for i in range(0, len(categoryList)):
  inputFileNameList = glob(inputLocation + categoryList[i] + "/*/" + prefix + "_tad.txt")
  for inputFileName in inputFileNameList: table[i].append(file_len(inputFileName) - 1)

# Writing table
maxN = max([len(e) for e in table])
outputFile = open(outputFileName, "w")
outputFile.write("\t".join(categoryList)+"\n")
for j in range(0, maxN):
  vec = []
  for i in range(0, len(table)):
    try: vec.append(str(table[i][j]))
    except Exception: vec.append("NA")
  outputFile.write("\t".join(vec)+"\n")
outputFile.close()


