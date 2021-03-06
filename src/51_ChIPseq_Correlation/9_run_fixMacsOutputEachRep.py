
#Import
import os
import sys
from glob import glob

# Macs list
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/8_Peaks_Separated_Replicate/"
tl = "/scratch/egadegu/FMOR/"
macsFileNameList = glob(il + "*_peaks.narrowPeak")

# Opening file
inFileName = fl + "9_fmor.txt"
inFile = open(inFileName, "w")

# Macs Loop
for macsFileName in macsFileNameList:

  # Parameter
  pkName = "_".join(macsFileName.split("_")[:-1])

  # Input
  percThreshold = "1"
  peakFileName = macsFileName
  summitFileName = pkName + "_summits.bed"
  tempLoc = tl + pkName + "/"
  outputPeakFileName = pkName + "_peaks_filtered.narrowPeak"
  outputSummitFileName = pkName + "_summits_filtered.bed"

  # Creating files
  inFile.write(" ".join([percThreshold, peakFileName, summitFileName, tempLoc, outputPeakFileName, outputSummitFileName])+"\n")

# Close
inFile.close()


