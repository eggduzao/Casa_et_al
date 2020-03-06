
# Import
import os
import sys

# Peak List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/34_Smc3_Analysis/input/"
pl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/34_Smc3_Analysis/0_Regions/"
bl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/7_FixedBigWig/"
tl = "/scratch/egadegu/CHS/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/34_Smc3_Analysis/1_Heatmap/"
peakList = ["Shared", "Stag1only", "Stag2only", "Rad21"]
peakHeight = ["48", "24", "3", "20"]
zMaxList = ["10", "5", "5", "5"]

# Open File
inputFileName = fl + "1_chs.txt"
inFile = open(inputFileName, "w")

# Peak Loop
for i in range(0, len(peakList)):

  # Bigwig List
  peakName = peakList[i]
  hei = peakHeight[i]
  bwList = ["Smc3"]

  # Bigwig Loop
  for bwName in bwList:

    # Input
    ext = "500"
    zMax = zMaxList[i]
    hheight = hei
    featureSummitFileName = pl + peakName + ".bed"
    bwFileName = bl + bwName + ".bw"
    bwLabel = peakName
    tempLocation = tl + peakName + "/"
    outputLocation = ol

    # Execution
    inFile.write(" ".join([ext, zMax, hheight, featureSummitFileName, bwFileName, bwLabel, tempLocation, outputLocation])+"\n")

inFile.close()


