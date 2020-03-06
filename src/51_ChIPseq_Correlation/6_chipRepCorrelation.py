
# Import
import os
import sys
from pysam import Samfile
from scipy.stats import spearmanr

###################################################################################################
# INPUT
###################################################################################################

# Input
regionFileName = sys.argv[1]
bamFileNameList = sys.argv[2].split(",")
outputFilePrefix = sys.argv[3]

###################################################################################################
# FUNCTIONS
###################################################################################################

def fetchTotalSignal(region, bamFile):
  res = 0
  try:
    for read in bamFile.fetch(region[0], int(region[1]), int(region[2])): res += 1
  except Exception: res = -1
  return res

###################################################################################################
# EXECUTION
###################################################################################################

# Create locations
outputLocation = "/".join(outputFilePrefix.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

# Open bam files
bamFileList = [Samfile(e, "rb") for e in bamFileNameList]

# Writing signal to table
vectorList = [[0.0] for e in range(0,len(bamFileList))]
regionFile = open(regionFileName, "rU")
for line in regionFile:
  ll = line.strip().split("\t")
  region = [ll[0], int(ll[1]), int(ll[2])]
  vector = []
  flag_success = True
  for bamFile in bamFileList:
    counts = fetchTotalSignal(region, bamFile)
    if(counts >= 0):
      vector.append(counts)
    else:
      flag_success = False
      break
  if(flag_success):
    for i in range(0,len(vector)): vectorList[i].append(vector[i])
regionFile.close()
for e in bamFileList: e.close()

# Calculating correlation
header = []
corr_matrix = [[[1.0, 1.0] for e in range(0,len(bamFileList))] for e in range(0,len(bamFileList))]
for i in range(0, len(bamFileList)-1):
  header.append(bamFileNameList[i].split("/")[-1].split(".")[0])
  for j in range(i+1, len(bamFileList)):
    rank, pvalue = spearmanr(vectorList[i], vectorList[j])
    corr_matrix[i][j] = [rank, pvalue]
for i in range(0, len(bamFileList)):
  for j in range(0, len(bamFileList)):
    if(i == j): corr_matrix[i][j] = [1.0, 1.0]
    elif(i > j): corr_matrix[i][j] = corr_matrix[j][i]
header.append(bamFileNameList[len(bamFileList)-1].split("/")[-1].split(".")[0])

# Writing matrix
outStatFileName = outputFilePrefix + "_stat.txt"
outPvalueFileName = outputFilePrefix + "_pvalue.txt"
outStatFile = open(outStatFileName, "w")
outPvalueFile = open(outPvalueFileName, "w")
outStatFile.write("\t".join(header)+"\n")
outPvalueFile.write("\t".join(header)+"\n")
for i in range(0, len(bamFileList)):
  outStatFile.write("\t".join([header[i]] + [str(e[0]) for e in corr_matrix[i]])+"\n")
  outPvalueFile.write("\t".join([header[i]] + [str(e[1]) for e in corr_matrix[i]])+"\n")
outStatFile.close()
outPvalueFile.close()


