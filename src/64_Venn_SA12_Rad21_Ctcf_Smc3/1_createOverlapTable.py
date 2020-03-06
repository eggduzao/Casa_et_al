
# Import
import os
import sys
import itertools

###################################################################################################
# INPUT
###################################################################################################

# Input
halfExt = int(sys.argv[1])
peakFileList = sys.argv[2].split(",")
temporaryLocation = sys.argv[3]
outputPairFileName = sys.argv[4]
outputTripleFileName = sys.argv[5]

###################################################################################################
# FUNCTIONS
###################################################################################################

def file_len(inFileName):
  i = -1
  with open(inFileName) as f:
    for i, l in enumerate(f): pass
  return i + 1

def extend(halfExt, inFileName, outFileName):
  inFile = open(inFileName, "rU")
  outFile = open(outFileName, "w")
  for line in inFile:
    ll = line.strip().split("\t")
    mid = (int(ll[1])+int(ll[2]))/2
    p1 = str(mid - halfExt); p2 = str(mid + halfExt)
    outFile.write("\t".join([ll[0], p1, p2] + ll[3:]) + "\n")
  inFile.close()
  outFile.close()

def sort(inFileName, outFileName):
  command = "sort -k1,1 -k2,2n "+inFileName+" > "+outFileName
  os.system(command)

def cut(inFileName, outFileName):
  command = "cut -f 1,2,3 "+inFileName+" > "+outFileName
  os.system(command)

def intersect(aFileName, bFileName, outFileName, intType = "u"):
  command = "intersectBed -"+intType+" -wa -a "+aFileName+" -b "+bFileName+" > "+outFileName
  os.system(command)

def pair_intersection(halfExt, inputFile1Name, inputFile2Name, temporaryLocation):

  # Cut file 1
  cutFileName1 = temporaryLocation + "cutFileName1.bed"
  cut(inputFile1Name, cutFileName1)

  # Cut file 2
  cutFileName2 = temporaryLocation + "cutFileName2.bed"
  cut(inputFile2Name, cutFileName2)

  # Extend file 1
  extFileName1 = temporaryLocation + "extFileName1.bed"
  extend(halfExt, cutFileName1, extFileName1)

  # Extend file 2
  extFileName2 = temporaryLocation + "extFileName2.bed"
  extend(halfExt, cutFileName2, extFileName2)

  # Sort file 1
  sortFileName1 = temporaryLocation + "sortFileName1.bed"
  sort(extFileName1, sortFileName1)

  # Sort file 2
  sortFileName2 = temporaryLocation + "sortFileName2.bed"
  sort(extFileName2, sortFileName2)

  # Intersection
  tempIntFileName = temporaryLocation + "tempIntFileName.bed"
  intersect(sortFileName1, sortFileName2, tempIntFileName, intType = "u")

  # area1, area2, n12
  return file_len(sortFileName1), file_len(sortFileName2), file_len(tempIntFileName)

def triple_intersection(halfExt, inputFile1Name, inputFile2Name, inputFile3Name, temporaryLocation):

  # Cut file 1
  cutFileName1 = temporaryLocation + "cutFileName1.bed"
  cut(inputFile1Name, cutFileName1)

  # Cut file 2
  cutFileName2 = temporaryLocation + "cutFileName2.bed"
  cut(inputFile2Name, cutFileName2)

  # Cut file 3
  cutFileName3 = temporaryLocation + "cutFileName3.bed"
  cut(inputFile3Name, cutFileName3)

  # Extend file 1
  extFileName1 = temporaryLocation + "extFileName1.bed"
  extend(halfExt, cutFileName1, extFileName1)

  # Extend file 2
  extFileName2 = temporaryLocation + "extFileName2.bed"
  extend(halfExt, cutFileName2, extFileName2)

  # Extend file 3
  extFileName3 = temporaryLocation + "extFileName3.bed"
  extend(halfExt, cutFileName3, extFileName3)

  # Sort file 1
  sortFileName1 = temporaryLocation + "sortFileName1.bed"
  sort(extFileName1, sortFileName1)

  # Sort file 2
  sortFileName2 = temporaryLocation + "sortFileName2.bed"
  sort(extFileName2, sortFileName2)

  # Sort file 3
  sortFileName3 = temporaryLocation + "sortFileName3.bed"
  sort(extFileName3, sortFileName3)

  # Intersection 1 2
  intFileName12 = temporaryLocation + "intFileName12.bed"
  intersect(sortFileName1, sortFileName2, intFileName12, intType = "u")

  # Intersection 2 3
  intFileName23 = temporaryLocation + "intFileName23.bed"
  intersect(sortFileName2, sortFileName3, intFileName23, intType = "u")

  # Intersection 1 3
  intFileName13 = temporaryLocation + "intFileName13.bed"
  intersect(sortFileName1, sortFileName3, intFileName13, intType = "u")

  # Intersection 1 2 3
  intFileName123 = temporaryLocation + "intFileName123.bed"
  intersect(intFileName12, sortFileName3, intFileName123, intType = "u")

  # area1, area2, area3, n12, n23, n13, n123
  return file_len(sortFileName1), file_len(sortFileName2), file_len(sortFileName3), file_len(intFileName12), file_len(intFileName23), file_len(intFileName13), file_len(intFileName123)

###################################################################################################
# EXECUTION
###################################################################################################

# Create output location
command = "mkdir -p "+temporaryLocation
os.system(command)
outputLocation = "/".join(outputPairFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

# Pairwise intersections
outputPairFile = open(outputPairFileName, "w")
for comb in list(itertools.combinations(peakFileList, 2)):
  inFileName1 = comb[0]
  inFileName2 = comb[1]
  inName1 = inFileName1.split("/")[-1].split(".")[0]
  inName2 = inFileName2.split("/")[-1].split(".")[0]
  area1, area2, n12 = pair_intersection(halfExt, inFileName1, inFileName2, temporaryLocation)
  outputPairFile.write("\t".join([inName1, inName2, str(area1), str(area2), str(n12)]) + "\n")
outputPairFile.close()

# Triple intersections
outputTripleFile = open(outputTripleFileName, "w")
for comb in list(itertools.combinations(peakFileList, 3)):
  inFileName1 = comb[0]
  inFileName2 = comb[1]
  inFileName3 = comb[2]
  inName1 = inFileName1.split("/")[-1].split(".")[0]
  inName2 = inFileName2.split("/")[-1].split(".")[0]
  inName3 = inFileName3.split("/")[-1].split(".")[0]
  area1, area2, area3, n12, n23, n13, n123 = triple_intersection(halfExt, inFileName1, inFileName2, inFileName3, temporaryLocation)
  outputTripleFile.write("\t".join([inName1, inName2, inName3, str(area1), str(area2), str(area3), str(n12), str(n23), str(n13), str(n123)]) + "\n")
outputTripleFile.close()


