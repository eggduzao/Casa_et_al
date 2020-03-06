
#Import
import os
import sys
from glob import glob

"""
# Macs list
counter = 1
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/input/"
il = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/4_Peaks/"
tl = "/scratch/egusmao/STAG_FMO/"
macsFileNameList = glob(il + "*_peaks.narrowPeak")

# Opening file
inFileName = fl + "5_fmo.txt"
inFile = open(inFileName, "w")

# Macs Loop
for macsFileName in macsFileNameList:

  # Input
  pValueThreshold = "9.0"
  peakFileName = macsFileName
  summitFileName = "_".join(macsFileName.split("_")[:-1]) + "_summits.bed"
  tempLoc = tl + str(counter) + "/"
  outputPeakFileName = "_".join(macsFileName.split("_")[:-1]) + "_peaks_filtered.narrowPeak"
  outputSummitFileName = "_".join(macsFileName.split("_")[:-1]) + "_summits_filtered.bed"

  # Creating files
  inFile.write(" ".join([pValueThreshold, peakFileName, summitFileName, tempLoc, outputPeakFileName, outputSummitFileName])+"\n")
  counter += 1

# Close
inFile.close()
"""



# Macs list
il = "/home/egg/Projects/Wendt_Stag/Data/tf/macs/"
tl = "./TEMP/"
macsFileNameList = [il+"HCT116_ChIP-seq_RAD21_HAIB_peaks.narrowPeak"]

# Macs Loop
for macsFileName in macsFileNameList:

  # Input
  pValueThreshold = "9.0"
  peakFileName = macsFileName
  summitFileName = "_".join(macsFileName.split("_")[:-1]) + "_summits.bed"
  tempLoc = tl
  outputPeakFileName = "_".join(macsFileName.split("_")[:-1]) + "_peaks_filtered.narrowPeak"
  outputSummitFileName = "_".join(macsFileName.split("_")[:-1]) + "_summits_filtered.bed"

  # Creating files
  command = "python 5_fixMacsOutput.py "+" ".join([pValueThreshold, peakFileName, summitFileName, tempLoc, outputPeakFileName, outputSummitFileName])
  os.system(command)


