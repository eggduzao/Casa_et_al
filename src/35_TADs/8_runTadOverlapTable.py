
# Import
import os
import sys

# Input
il = "/home/egg/Projects/Wendt_Stag/Results/4_TADs/7_1_Individual_TADs_Bed_Files/"
resolution = "25000"
allowedShiftBins = "0"
tadLabelList = ",".join(["STAG1_WT", "STAG1_DEG", "STAG2_WT", "STAG2_DEG"])
tadFileNameList = ",".join([il + e + ".bed" for e in ["STAG1_WT", "STAG1_DEG", "STAG2_WT", "STAG2_DEG"]])
tempLoc = "./TEMP/"
outputFileName = "/home/egg/Projects/Wendt_Stag/Results/4_TADs/8_Tad_Overlap/overlap.txt"

# Execution
command = "python 8_tadOverlapTable.py "+" ".join([resolution, allowedShiftBins, tadLabelList, tadFileNameList, tempLoc, outputFileName])
os.system(command)


