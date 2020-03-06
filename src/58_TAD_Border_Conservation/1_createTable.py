
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from glob import glob

# Input
resolution = int(sys.argv[1])
threshold = int(sys.argv[2])
tadLocation = sys.argv[3]
tadNamePlusList = sys.argv[4].split(",")
tadNameMinusList = sys.argv[5].split(",")
outputFileName = sys.argv[6]

# Initialization
outputLocation = "/".join(outputFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)
tadFileNamePlusList = [tadLocation + e + "_tad.txt" for e in tadNamePlusList]
tadFileNameMinusList = [tadLocation + e + "_tad.txt" for e in tadNameMinusList]

###################################################################################################
# Functions
###################################################################################################

def round_down(num, divisor):
  return num - (num%divisor)

def round_up(num, divisor):
  return num + (divisor-(num%divisor))

def create_tad_dict(threshold, resolution, tadFileName):

  # Create TAD dict
  tad_dict = dict()
  tadFile = open(tadFileName, "rU")
  tadFile.readline()
  for line in tadFile:
    ll = line.strip().split("\t")
    t1 = round_down(int(ll[0]), 25000)
    for i in range(t1 - (threshold * resolution), t1 + (threshold * resolution), resolution):
      tad_dict[i] = True
  tadFile.close()

  # Return objects
  return tad_dict

###################################################################################################
# Execution
###################################################################################################

outputFile = open(outputFileName, "w")
for i in range(0, len(tadFileNamePlusList)):
  tadFileNamePlus = tadFileNamePlusList[i]
  tadFileNameMinus = tadFileNameMinusList[i]
  tad_dict = create_tad_dict(threshold, resolution, tadFileNamePlus)
  tadFile = open(tadFileNameMinus, "rU")
  tadFile.readline()
  summ = 0.0
  count = 0.0
  for line in tadFile:
    ll = line.strip().split("\t")
    try:
      a = tad_dict[round_down(int(ll[0]), 25000)]
      summ += 1.0
    except Exception: pass
    count += 1.0
  tadFile.close()
  outputFile.write(str(round(summ/count,2))+"\n")
outputFile.close()



