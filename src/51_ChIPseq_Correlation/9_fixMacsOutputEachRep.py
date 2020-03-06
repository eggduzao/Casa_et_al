
#Import
import os
import sys
from scipy.stats import percentileofscore

###################################################################################################
# Input
###################################################################################################

# Input
percThreshold = int(sys.argv[1])
peakFileName = sys.argv[2]
summitFileName = sys.argv[3]
tempLoc = sys.argv[4]
outputPeakFileName = sys.argv[5]
outputSummitFileName = sys.argv[6]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

###################################################################################################
# Execution
###################################################################################################

# Allowed chromosomes
chrDict = dict([("chr"+str(e), True) for e in range(1,23)+["X"]])

# 1. Creating dictionary of filtered (by p-value) peak names
# 2. Creating list of p-values.
pValueList = []
peakDict = dict() # PEAK_NAME -> [ [<a>, all the other narrowpeak fields in order], [<b>, all the other narrowpeak fields in order], ...]
peakFile = open(peakFileName, "rU")
for line in peakFile:
  ll = line.strip().split("\t")
  try: valid_chromosome = chrDict[ll[0]]
  except Exception: continue
  pValue = float(ll[7])
  lastName = ll[3].split("_")[-1]
  lastChar = lastName[-1]
  try:
    lastChar = int(lastChar)
    key = ll[3]
    ll = ["."] + ll
  except Exception:
    key = ll[3][:-1]
    ll = [ll[3][-1]] + ll
  try: peakDict[key].append(ll)
  except Exception: peakDict[key] = [ll]
  pValueList.append(pValue)
peakFile.close()

# 1. Writing only the best peaks for each replicate (a, b, c, etc...)
# 2. Creating dictionary to filter summits
summitDict = dict() # PEAK_NAME -> PERCENTILE
tempPeakFileName = tempLoc + "tempPeakFileName.bed"
outputPeakFile = open(tempPeakFileName, "w")
for key in peakDict.keys():
  peakList = peakDict[key]
  bestSuffix = None
  bestPeak = None
  bestScore = -1
  for peak in peakList:
    currScore = float(peak[8])
    if(currScore > bestScore):
      if(peak[0] == "."): bestSuffix = ""
      else: bestSuffix = peak[0]
      bestPeak = peak[1:]
      bestScore = currScore
  perc = percentileofscore(pValueList, bestScore)
  percentile = str(perc)
  if(perc < percThreshold): continue
  outputPeakFile.write("\t".join(bestPeak[:3] + [key] + bestPeak[4:] + [percentile])+"\n")
  summitDict[key+bestSuffix] = percentile
outputPeakFile.close()

# Writing only allowed summits
summitFile = open(summitFileName, "rU")
tempSummitFileName = tempLoc + "tempSummitFileName.bed"
outputSummitFile = open(tempSummitFileName, "w")
for line in summitFile:
  ll = line.strip().split("\t")
  try: percentile = summitDict[ll[3]]
  except Exception: continue
  lastName = ll[3].split("_")[-1]
  lastChar = lastName[-1]
  try:
    lastChar = int(lastChar)
    key = ll[3]
  except Exception: key = ll[3][:-1]
  outputSummitFile.write("\t".join(ll[:3] + [key, percentile, "."])+"\n")
summitFile.close()
outputSummitFile.close()

# Sorting peak file
command = "sort -k1,1 -k2,2n "+tempPeakFileName+" > "+outputPeakFileName
os.system(command)

# Sorting summit file
command = "sort -k1,1 -k2,2n "+tempSummitFileName+" > "+outputSummitFileName
os.system(command)

# Initialization
command = "rm -rf "+tempLoc
os.system(command)


