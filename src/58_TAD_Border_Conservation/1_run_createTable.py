
# Import
import os
import sys

# TAD List
prefix = "T_3_80_25_100_5_10_0.95_0.5"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/29_TAD_Border_Conservation/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/4_TADs/1_TADs_GMAP/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/29_TAD_Border_Conservation/"
tadList = [["69_127-", "69_127plus"], ["3B9_5-", "3B9_5plus"]]
nameList = ["STAG1", "STAG2"]

# Open File
inputFileName = fl + "1_ctc.txt"
inFile = open(inputFileName, "w")

chrList = ["chr" + str(e) for e in range(1,23) + ["X"]]

# TAD Loop
for i in range(0,len(tadList)):

  # Parameter
  name = nameList[i]
  plusFile = tadList[i][1]
  minusFile = tadList[i][0]

  # Input
  resolution = "25000"
  threshold = "2"
  tadLocation = il
  tadNamePlusList = ",".join([plusFile + "/" + e + "/" + prefix for e in chrList])
  tadNameMinusList = ",".join([minusFile + "/" + e + "/" + prefix for e in chrList])
  outputFileName = ol + name + ".txt"

  # Execution
  inFile.write(" ".join([resolution, threshold, tadLocation, tadNamePlusList, tadNameMinusList, outputFileName])+"\n")

inFile.close()


