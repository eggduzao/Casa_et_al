
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

# Input
halfExt = int(sys.argv[1])
totalCount = float(sys.argv[2])
regionFileName = sys.argv[3]
bamFileName = sys.argv[4]
outputFileName = sys.argv[5]

# Initialization
outputLocation = "/".join(outputFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

# Fetch total signal
def fetchTotalSignal(valueToAdd, region, bamFile):
  res = 0.0
  try:
    for read in bamFile.fetch(region[0], int(region[1]), int(region[2])): res += valueToAdd
  except Exception: pass
  return res

###################################################################################################
# Execution
###################################################################################################

# Values
nValue = 1.0
rpmValue = 1000000.0 / totalCount

# Iterating in regions
resTable = []
regionFile = open(regionFileName)
bamFile = Samfile(bamFileName)
outputFile = open(outputFileName, "w")
for line in regionFile:
  ll = line.strip().split("\t")
  mid = (int(ll[1]) + int(ll[2])) / 2
  p1 = str(mid - halfExt); p2 = str(mid + halfExt)
  try: signalN = fetchTotalSignal(nValue, [ll[0], p1, p2], bamFile)
  except Exception: signalN = 0.0
  try: signalRpm = fetchTotalSignal(rpmValue, [ll[0], p1, p2], bamFile)
  except Exception: signalRpm = 0.0
  outputFile.write("\t".join([ll[0], p1, p2, str(signalN), str(signalRpm)])+"\n")
regionFile.close()
bamFile.close()
outputFile.close()


