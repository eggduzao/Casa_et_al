
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
abTreatFileName = sys.argv[1]
abContrFileName = sys.argv[2]
outputFileName = sys.argv[3]
outputBedFilePrefix = sys.argv[4]

# Creating treatment dictionary
abTreatFile = open(abTreatFileName, "rU")
abTreatDict = dict()
for line in abTreatFile:
  ll = line.strip().split("\t")
  key = ":".join(ll[:3])
  abTreatDict[key] = float(ll[4])
abTreatFile.close()

# Creating table
abContrFile = open(abContrFileName, "rU")
outputFile = open(outputFileName, "w")
outputFile.write("\t".join(["X", "Y", "Z"])+"\n")
for line in abContrFile:
  ll = line.strip().split("\t")
  key = ":".join(ll[:3])
  try: treatValue = abTreatDict[key]
  except Exception: continue
  contrValue = float(ll[4])
  if(treatValue < 0 and contrValue < 0): sign = "1"
  elif(treatValue > 0 and contrValue < 0): sign = "2"
  elif(treatValue < 0 and contrValue > 0): sign = "3"
  else: sign = "4"
  outputFile.write("\t".join([str(contrValue), str(treatValue), sign])+"\n")
abContrFile.close()
outputFile.close()

# Creating bed files
abContrFile = open(abContrFileName, "rU")
aaBedFileName = outputBedFilePrefix+"_AA.bed"
aaBedFile = open(aaBedFileName, "w")
abBedFileName = outputBedFilePrefix+"_AB.bed"
abBedFile = open(abBedFileName, "w")
baBedFileName = outputBedFilePrefix+"_BA.bed"
baBedFile = open(baBedFileName, "w")
bbBedFileName = outputBedFilePrefix+"_BB.bed"
bbBedFile = open(bbBedFileName, "w")
for line in abContrFile:
  ll = line.strip().split("\t")
  key = ":".join(ll[:3])
  try: treatValue = abTreatDict[key]
  except Exception: continue
  contrValue = float(ll[4])
  if(treatValue < 0 and contrValue < 0):
    outFile = bbBedFile
    name = "BB"
  elif(treatValue > 0 and contrValue < 0):
    outFile = baBedFile
    name = "BA"
  elif(treatValue < 0 and contrValue > 0):
    outFile = abBedFile
    name = "AB"
  else:
    outFile = aaBedFile
    name = "AA"
  outFile.write("\t".join(ll[:3]+[name])+"\n")
abContrFile.close()
aaBedFile.close()
baBedFile.close()
abBedFile.close()
bbBedFile.close()


