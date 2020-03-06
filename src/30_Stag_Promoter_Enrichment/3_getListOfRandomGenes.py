
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
inputActiveFileName = sys.argv[1]
inputInactiveFileName = sys.argv[2]
activeGenesFileName = sys.argv[3]
inactiveGenesFileName = sys.argv[4]
tempLocation = sys.argv[5]
outputActiveControlFileName = sys.argv[6]
outputInactiveControlFileName = sys.argv[7]

# Initialization
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def file_len(fname):
  i = -1
  with open(fname) as f:
    for i, l in enumerate(f): pass
  return i + 1
  
###################################################################################################
# Execution
###################################################################################################

# Number of genes is 100 times the number of genes in the treatment list
numberOfGenesA = 100 * file_len(inputActiveFileName)
numberOfGenesI = 100 * file_len(inputInactiveFileName)

# Get random genes
tempFileAName = tempLocation + "tempFileAName.bed"
command = "shuf -n "+str(numberOfGenesA)+" "+activeGenesFileName+" > "+tempFileAName
os.system(command)
tempFileIName = tempLocation + "tempFileIName.bed"
command = "shuf -n "+str(numberOfGenesI)+" "+inactiveGenesFileName+" > "+tempFileIName
os.system(command)

# Sorting lists
command = "sort -k1,1 -k2,2n "+tempFileAName+" > "+outputActiveControlFileName
os.system(command)
command = "sort -k1,1 -k2,2n "+tempFileIName+" > "+outputInactiveControlFileName
os.system(command)

# Removing temp
command = "rm -rf "+tempLocation
os.system(command)


