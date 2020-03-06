
# Import
import os
import sys

# Input
chrom = sys.argv[1]
resolution = int(sys.argv[2])
chromSizesFileName = sys.argv[3]
matrixFileName = sys.argv[4]
outputFileName = sys.argv[5]

# Initialization
outputLocation = "/".join(matrixFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)
outputLocation = "/".join(outputFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def read_chrom_sizes(chromSizesFileName):

  chromSizesDict = dict()
  chromSizesFile = open(chromSizesFileName,"rU")
  for line in chromSizesFile:
    ll = line.strip().split("\t")
    chromSizesDict[ll[0]] = int(ll[1])
  chromSizesFile.close()

  return chromSizesDict

def read_matrix_dictionary(chrom, matrixFileName):

  matrixDict = dict()
  matrixFile = open(matrixFileName, "rU")
  for line in matrixFile:
    ll = line.strip().split("\t")
    if(ll[0] != chrom): continue
    key1 = ":".join([ll[0],ll[1],ll[2]]); key2 = ":".join([ll[0],ll[2],ll[1]])
    matrixDict[key1] = ll[3]
    matrixDict[key2] = ll[3]
  matrixFile.close()

  return matrixDict

def create_full_matrix(chrom, resolution, chromSizesFileName, matrixFileName, outputFileName):

  # Read chrom sizes
  chromSizesDict = read_chrom_sizes(chromSizesFileName)

  # Read matrix dictionary
  matrixDict = read_matrix_dictionary(chrom, matrixFileName)

  # Creating the full matrix
  outputFile = open(outputFileName, "w")
  for i in range(0, chromSizesDict[chrom] + resolution, resolution):
    if(i >= chromSizesDict[chrom]): continue
    vec = []
    for j in range(0, chromSizesDict[chrom] + resolution, resolution):
      if(j >= chromSizesDict[chrom]): continue
      key = ":".join([chrom, str(i), str(j)])
      value = "0"
      try: value = matrixDict[key]
      except Exception: pass
      vec.append(value)
    outputFile.write("\t".join(vec)+"\n")
  outputFile.close()      

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
create_full_matrix(chrom, resolution, chromSizesFileName, matrixFileName, outputFileName)


