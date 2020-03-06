
#Import
import os
import sys

###################################################################################################
# INPUT
###################################################################################################

# Input
intType = int(sys.argv[1]) # 1 (only one is enough) or 2 (has to overlap both) sides of the loop
resolution = int(sys.argv[2])
halfExt = int(sys.argv[3])
inputLoopNameList = sys.argv[4].split(",")
inputLoopFileNameList = sys.argv[5].split(",")
inputPeaksNameList = sys.argv[6].split(",")
inputPeaksFileNameList = sys.argv[7].split(",")
temporaryLocation = sys.argv[8]
outputFileName = sys.argv[9]

# Create temp and output location
command = "mkdir -p "+temporaryLocation
os.system(command)
outputLocation = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# FUNCTIONS
###################################################################################################

def cut(inFileName, outFileName):
  command = "cut -f 1,2,3 "+inFileName+" > "+outFileName
  os.system(command)

def extend_region(halfExt, inFileName, outFileName):
  inFile = open(inFileName, "rU")
  outFile = open(outFileName, "w")
  for line in inFile:
    ll = line.strip().split("\t")
    mid = (int(ll[1])+int(ll[2]))/2
    p1 = str(mid - halfExt); p2 = str(mid + halfExt)
    outFile.write("\t".join([ll[0], p1, p2] + ll[3:]) + "\n")
  inFile.close()
  outFile.close()

def extend_loop(resolution, inFileName, outFileName):
  inFile = open(inFileName, "rU")
  outFile = open(outFileName, "w")
  for line in inFile:
    ll = line.strip().split("\t")
    a1 = str(int(ll[1])); a2 = str(int(ll[1]) + resolution)
    b1 = str(int(ll[2])); b2 = str(int(ll[2]) + resolution)
    outFile.write("\t".join([ll[0], a1, a2, b1, b2]) + "\n")
  inFile.close()
  outFile.close()

def intersect(lFileName, rFileName, outFileName, intType = 1):
  rDict = dict()
  rFile = open(rFileName, "rU")
  for line in rFile:
    ll = line.strip().split("\t")
    for i in range(int(ll[1]), int(ll[2])): rDict[i] = True
  rFile.close()
  lFile = open(lFileName, "rU")
  outFile = open(outFileName, "w")
  for line in lFile:
    ll = line.strip().split("\t")
    ovl1 = False; ovl2 = False
    for i in range(int(ll[1]), int(ll[2])):
      try:
        ovl1 = rDict[i]
        break
      except Exception: pass
    for i in range(int(ll[3]), int(ll[4])):
      try:
        ovl2 = rDict[i]
        break
      except Exception: pass
    if(intType == 1 and (ovl1 or ovl2)): outFile.write("\t".join([ll[0], ll[1], ll[3]])+"\n")
    elif(intType == 2 and (ovl1 and ovl2)): outFile.write("\t".join([ll[0], ll[1], ll[3]])+"\n")
  rDict = None
  lFile.close()
  outFile.close()

def get_distribution(intType, resolution, halfExt, loopFileName, regionsFileName, temporaryLocation):

  # Cut files
  loopFileNameCut = temporaryLocation + "loopFileNameCut.bed"
  cut(loopFileName, loopFileNameCut)
  regionsFileNameCut = temporaryLocation + "regionsFileNameCut.bed"
  cut(regionsFileName, regionsFileNameCut)

  # Extend files
  loopFileNameExt = temporaryLocation + "loopFileNameExt.bed"
  extend_loop(resolution, loopFileNameCut, loopFileNameExt)
  regionsFileNameExt = temporaryLocation + "regionsFileNameExt.bed"
  extend_region(halfExt, regionsFileNameCut, regionsFileNameExt)

  # Intersect loops and regions
  intFileName = temporaryLocation + "intFileName.bed"
  intersect(loopFileNameExt, regionsFileNameExt, intFileName, intType = intType)

  # Getting distribution
  resVec = []
  intFile = open(intFileName, "rU")
  for line in intFile:
    ll = line.strip().split("\t")
    resVec.append(str(int(ll[2]) - int(ll[1])))
  intFile.close()  

  return resVec

###################################################################################################
# EXECUTION
###################################################################################################

# Iterating on input files to fetch distribution
distMat = []
header = []
for i in range(0, len(inputLoopFileNameList)):
  for j in range(0, len(inputPeaksFileNameList)):
    header.append("_".join([inputLoopNameList[i], inputPeaksNameList[j]]))
    distMat.append(get_distribution(intType, resolution, halfExt, inputLoopFileNameList[i], inputPeaksFileNameList[j], temporaryLocation))

# Writing output
maxL = max([len(e) for e in distMat])
outputFile = open(outputFileName, "w")
outputFile.write("\t".join(header)+"\n")
for j in range(0, maxL):
  vec = []
  for i in range(0, len(distMat)):
    try: vec.append(distMat[i][j])
    except Exception: vec.append("NA")
  outputFile.write("\t".join(vec)+"\n")
outputFile.close()


