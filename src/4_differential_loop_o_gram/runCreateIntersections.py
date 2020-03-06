
# Import
import os
import sys

# Input files
ml = "/home/egg/Projects/Papantonis_Stag/Results/loop_o_gram/matrix/"
plusFileList = [ml+"3B9_5plus:STAG1only_regions.txt", ml+"69_127plus:STAG1only_regions.txt"]
minusFileList = [ml+"3B9_5-:STAG1only_regions.txt", ml+"69_127-:STAG1only_regions.txt"]
outLoc = "/home/egg/Projects/Papantonis_Stag/Results/differential_loop_o_gram/intersections/"

# Loop on hic files
for i in range(0,len(plusFileList)):

  # Parameters
  stagName = "STAG1"
  if("STAG2" in plusFileList[i]): stagName = "STAG2"
  expName = "127"
  if("3B9" in plusFileList[i]): expName = "3B9"

  # Input
  plusRegionFileName = plusFileList[i]
  minusRegionFileName = minusFileList[i]
  tempLocation = outLoc+"TEMP/"
  plusOnlyOutputFileName = outLoc+"_".join([stagName, expName, "plusOnly"])+".bed"
  minusOnlyOutputFileName = outLoc+"_".join([stagName, expName, "minusOnly"])+".bed"
  intersectionOutputFileName = outLoc+"_".join([stagName, expName, "intersection"])+".bed"

  # Execution
  command = "python createIntersections.py "+" ".join([plusRegionFileName, minusRegionFileName, tempLocation, plusOnlyOutputFileName, minusOnlyOutputFileName, intersectionOutputFileName])
  os.system(command)


