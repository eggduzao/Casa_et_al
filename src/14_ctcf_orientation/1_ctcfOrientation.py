
# Import
import os
import sys

###################################################################################################
# Input
###################################################################################################

# Input
loopFileName = sys.argv[1]
motifFileName = sys.argv[2]
tempLocation = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def getBestMotif(region, motifFileName, tempLocation):
  regionFileName = tempLocation+"regionFileName.bed"
  regionFile = open(regionFileName, "w")
  regionFile.write("\t".join(region)+"\n")
  regionFile.close()
  intFileName = tempLocation+"intFileName.bed"
  command = "intersectBed -wa -u -a "+motifFileName+" -b "+regionFileName+" > "+intFileName
  os.system(command)
  intFile = open(intFileName, "rU")
  currMotif = None
  currScore = -999999.9
  for line in intFile:
    ll = line.strip().split("\t")
    if(float(ll[4]) > currScore):
      currScore = float(ll[4])
      currMotif = ll[1:]
  intFile.close()
  return currMotif
  

###################################################################################################
# Execution
###################################################################################################

# Iterating on loop anchor's file
loopFile = open(loopFileName, "rU")
outputFile = open(outputFileName, "w")
for line in loopFile:

  # Initialization
  ll = line.strip().split("\t")
  chrom = ll[0]; x1 = ll[1]; x2 = ll[2]; y1 = ll[4]; y2 = ll[5]
  
  # Fetching highest scored motifs within loop anchors
  m1 = getBestMotif([chrom, x1, x2], motifFileName, tempLocation)
  m2 = getBestMotif([chrom, y1, y2], motifFileName, tempLocation)

  # Creating and writing vector
  vector = [chrom, x1, x2]
  if(m1):
    if(m1[4] == "+"): m1[4] = "P"
    else: m1[4] = "N"
    vector = vector + [m1[0], m1[1], m1[2], m1[3], m1[4]]
  else: vector = vector + ["NA", "NA", "NA","NA", "NA"]
  vector = vector + [y1, y2]
  if(m2):
    if(m2[4] == "+"): m2[4] = "P"
    else: m2[4] = "N"
    vector = vector + [m2[0], m2[1], m2[2], m2[3], m2[4]]
  else: vector = vector + ["NA", "NA", "NA","NA", "NA"]
  outputFile.write("\t".join(vector)+"\n")

# Closing files
loopFile.close()
outputFile.close()

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


