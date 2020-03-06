
# Import
import os
import sys

# Input
chrom = sys.argv[1]
pos1 = int(sys.argv[2])
pos2 = int(sys.argv[3])
resolution = int(sys.argv[4])
matrixFileName = sys.argv[5]
outputFileName = sys.argv[6]

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

def read_matrix_dictionary(chrom, pos1, pos2, matrixFileName):

  matrixDict = dict()
  matrixFile = open(matrixFileName, "rU")
  for line in matrixFile:
    ll = line.strip().split("\t")
    if(ll[0] != chrom or pos1 > int(ll[1]) or pos2 < int(ll[2])): continue
    key1 = ":".join([ll[0],ll[1],ll[2]]); key2 = ":".join([ll[0],ll[2],ll[1]])
    matrixDict[key1] = ll[3]
    matrixDict[key2] = ll[3]
  matrixFile.close()

  return matrixDict

def create_full_matrix(chrom, pos1, pos2, resolution, matrixFileName, outputFileName):

  # Read matrix dictionary
  matrixDict = read_matrix_dictionary(chrom, pos1, pos2, matrixFileName)

  # Creating the full matrix
  outputFile = open(outputFileName, "w")
  for i in range(pos1, pos2 + resolution, resolution):
    if(i >= pos2): continue
    vec = []
    for j in range(pos1, pos2 + resolution, resolution):
      if(j >= pos2): continue
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
create_full_matrix(chrom, pos1, pos2, resolution, matrixFileName, outputFileName)


