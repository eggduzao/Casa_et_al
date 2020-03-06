
# Import
import os
import sys

# Resolution List
ml = "/home/egg/Projects/Wendt_Stag/Previous_Results/3_Loopograms/matrix_new/"
outLoc = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/input/"
resList = ["10K_none", "10K_norm", "25K_none", "25K_norm", "50K_none", "50K_norm"]
resList = ["25K_norm"]

# Resolution Loop
for res in resList:

  # Input Lists
  plusFileList = [ml+res+"/69_127plus.txt:STAG1_regions.txt", ml+res+"/3B9_5plus.txt:STAG2_regions.txt"]
  minusFileList = [ml+res+"/69_127-.txt:STAG1_regions.txt", ml+res+"/3B9_5-.txt:STAG2_regions.txt"]

  # Input Loop
  for i in range(0,len(plusFileList)):

    # Parameters
    stagName = "R_STAG1"
    if("STAG2" in plusFileList[i]): stagName = "R_STAG2"
    expName = "RM_STAG1"
    if("3B9" in plusFileList[i]): expName = "RM_STAG2"

    # Input
    plusRegionFileName = plusFileList[i]
    minusRegionFileName = minusFileList[i]
    tempLocation = outLoc+res+"/TEMP/"
    plusOnlyOutputFileName = outLoc+res+"/"+"_".join([stagName, expName, "plusOnly"])+".bed"
    minusOnlyOutputFileName = outLoc+res+"/"+"_".join([stagName, expName, "minusOnly"])+".bed"
    intersectionOutputFileName = outLoc+res+"/"+"_".join([stagName, expName, "intersection"])+".bed"

    # Execution
    command = "python 0_createIntersectionsAndUnion.py "+" ".join([plusRegionFileName, minusRegionFileName, tempLocation, plusOnlyOutputFileName, minusOnlyOutputFileName, intersectionOutputFileName])
    os.system(command)


