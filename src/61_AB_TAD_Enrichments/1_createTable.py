
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from math import log

# Input
resolution = int(sys.argv[1])
regionFileName = sys.argv[2]
abContrFileName = sys.argv[3]
abTreatFileName = sys.argv[4]
tadContrFileName = sys.argv[5]
tadTreatFileName = sys.argv[6]
matrixContrFileName = sys.argv[7]
matrixTreatFileName = sys.argv[8]
temporaryLocation = sys.argv[9]
outputFileName = sys.argv[10]

# Initialization
command = "mkdir -p "+temporaryLocation
os.system(command)
outputLocation = "/".join(outputFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def round_down(num, divisor):
  return num - (num%divisor)

def round_up(num, divisor):
  return num + (num%divisor)

def read_matrix_dictionary(matrix_file_name):

  # Initialization
  matrix_dict = dict()
  matrix_file = open(matrix_file_name,"rU")
  bad_list = ["NaN", "Inf", "-Inf", "NA", "0", "0.0"]

  # Populating matrix
  for line in matrix_file:
    ll = line.strip().split("\t")
    if(ll[3] in bad_list): continue
    p1 = int(ll[1])
    p2 = int(ll[2])
    minP = str(min(p1, p2))
    maxP = str(max(p1, p2))
    matrix_dict[":".join([ll[0], minP, maxP])] = float(ll[3])
    matrix_dict[":".join([ll[0], maxP, minP])] = float(ll[3])

  # Returning objects
  matrix_file.close()
  return matrix_dict

def fetchTotalSignal(resolution, matrixDict, regionFileName):
  total = 0.0
  regionFile = open(regionFileName, "rU")
  for line in regionFile:
    ll = line.strip().split()
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])
    for i in range(round_down(p1, resolution), round_up(p2, resolution), resolution):
      key = ":".join([chrom, str(i), str(i + resolution)])
      try: total += matrixDict[key]
      except Exception: pass
  regionFile.close()
  return total

def intersection(fileName1, fileName2, intFileName):
  command = "intersectBed -u -wa -a "+fileName1+" -b "+fileName2+" > "+intFileName
  os.system(command)

def splitAB(abFileName, aFileName, bFileName, abBorderFileName):

  abFile = open(abFileName, "rU")
  aFile = open(aFileName, "w")
  bFile = open(bFileName, "w")
  abBorderFile = open(abBorderFileName, "w")
  prevChrom = "X"
  prevStart = 0
  prevComp = "X"
  for line in abFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = ll[1]; p2 = ll[2]; comp = ll[3]
    if(chrom != prevChrom):
      prevChrom = chrom
      prevStart = p1
      prevComp = comp
      continue
    if(comp != prevComp):
      if(prevComp == "A"): aFile.write("\t".join([chrom, prevStart, p1])+"\n")
      elif(prevComp == "B"): bFile.write("\t".join([chrom, prevStart, p1])+"\n")
      abBorderFile.write("\t".join([chrom, str(int(p1) - 1000), str(int(p1) + 1000)])+"\n")
      prevStart = p1
      prevComp = comp
  abFile.close()
  aFile.close()
  bFile.close()
  abBorderFile.close()

def getTadBorders(tadFileName, tadBorderFileName):

  tadFile = open(tadFileName, "rU")
  tadBorderFile = open(tadBorderFileName, "w")
  prevTad = tadFile.readline().strip().split("\t")
  prevChrom = prevTad[0]
  for line in tadFile:
    ll = line.strip().split("\t")
    if(ll[0] != prevChrom):
      prevTad = ll
      prevChrom = ll[0]
      continue
    mid = int(ll[1])
    tadBorderFile.write("\t".join([ll[0], str(mid - 1000), str(mid + 1000)])+"\n")
  tadFile.close()
  tadBorderFile.close()

###################################################################################################
# Execution
###################################################################################################

# Read matrices
matrixContr = read_matrix_dictionary(matrixContrFileName)
matrixTreat = read_matrix_dictionary(matrixTreatFileName)

# AB Files
aContrFileName = temporaryLocation + "aContrFileName.bed"
bContrFileName = temporaryLocation + "bContrFileName.bed"
abBorderContrFileName = temporaryLocation + "abBorderContrFileName.bed"
splitAB(abContrFileName, aContrFileName, bContrFileName, abBorderContrFileName)
aTreatFileName = temporaryLocation + "aTreatFileName.bed"
bTreatFileName = temporaryLocation + "bTreatFileName.bed"
abBorderTreatFileName = temporaryLocation + "abBorderTreatFileName.bed"
splitAB(abTreatFileName, aTreatFileName, bTreatFileName, abBorderTreatFileName)

# TAD Files
tadBorderContrFileName = temporaryLocation + "tadBorderContrFileName.bed"
getTadBorders(tadContrFileName, tadBorderContrFileName)
tadBorderTreatFileName = temporaryLocation + "tadBorderTreatFileName.bed"
getTadBorders(tadTreatFileName, tadBorderTreatFileName)

# Calculating intersections
regInAContrFileName = temporaryLocation + "regInAContrFileName.bed"
intersection(regionFileName, aContrFileName, regInAContrFileName)
regInATreatFileName = temporaryLocation + "regInATreatFileName.bed"
intersection(regionFileName, aTreatFileName, regInATreatFileName)
regInBContrFileName = temporaryLocation + "regInBContrFileName.bed"
intersection(regionFileName, bContrFileName, regInBContrFileName)
regInBTreatFileName = temporaryLocation + "regInBTreatFileName.bed"
intersection(regionFileName, bTreatFileName, regInBTreatFileName)
regInABborderContrFileName = temporaryLocation + "regInABborderContrFileName.bed"
intersection(regionFileName, abBorderContrFileName, regInABborderContrFileName)
regInABborderTreatFileName = temporaryLocation + "regInABborderTreatFileName.bed"
intersection(regionFileName, abBorderTreatFileName, regInABborderTreatFileName)
regInTADContrFileName = temporaryLocation + "regInTADContrFileName.bed"
intersection(regionFileName, tadBorderContrFileName, regInTADContrFileName)
regInTADTreatFileName = temporaryLocation + "regInTADTreatFileName.bed"
intersection(regionFileName, tadBorderTreatFileName, regInTADTreatFileName)

# Writing output
outputFile = open(outputFileName, "w")
outputFile.write("A\t" + str(round((fetchTotalSignal(resolution, matrixTreat, regInATreatFileName)+1) / (fetchTotalSignal(resolution, matrixContr, regInAContrFileName)+1),4)) + "\t")
outputFile.write("B\t" + str(round((fetchTotalSignal(resolution, matrixTreat, regInBTreatFileName)+1) / (fetchTotalSignal(resolution, matrixContr, regInBContrFileName)+1),4)) + "\t")
outputFile.write("AB\t" + str(round((fetchTotalSignal(resolution, matrixTreat, regInABborderTreatFileName)+1) / (fetchTotalSignal(resolution, matrixContr, regInABborderContrFileName)+1),4)) + "\t")
outputFile.write("TAD\t" + str(round((fetchTotalSignal(resolution, matrixTreat, regInTADTreatFileName)+1) / (fetchTotalSignal(resolution, matrixContr, regInTADContrFileName)+1),4)) + "\n")
outputFile.write("A\t" + str(round(log(((fetchTotalSignal(resolution, matrixTreat, regInATreatFileName)+1) / (fetchTotalSignal(resolution, matrixContr, regInAContrFileName)+1))+1),4)) + "\t")
outputFile.write("B\t" + str(round(log(((fetchTotalSignal(resolution, matrixTreat, regInBTreatFileName)+1) / (fetchTotalSignal(resolution, matrixContr, regInBContrFileName)+1))+1),4)) + "\t")
outputFile.write("AB\t" + str(round(log(((fetchTotalSignal(resolution, matrixTreat, regInABborderTreatFileName)+1) / (fetchTotalSignal(resolution, matrixContr, regInABborderContrFileName)+1))+1),4)) + "\t")
outputFile.write("TAD\t" + str(round(log(((fetchTotalSignal(resolution, matrixTreat, regInTADTreatFileName)+1) / (fetchTotalSignal(resolution, matrixContr, regInTADContrFileName)+1))+1),4)) + "\t")
outputFile.close()


