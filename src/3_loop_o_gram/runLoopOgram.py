
# Import
import os
import sys

# Location List
counter = 1
il = "/home/egg/Projects/Wendt_Stag/Data/stag_matrix_files/"
rl = "/home/egg/Projects/Wendt_Stag/Data/stag_bed_files/"
olm = "/home/egg/Projects/Wendt_Stag/Previous_Results/3_Loopograms/matrix_new/"
olg = "/home/egg/Projects/Wendt_Stag/Previous_Results/3_Loopograms/graph_new/"
#locationList = ["10K_none", "10K_norm", "25K_none", "25K_norm", "50K_none", "50K_norm"]

locationList = ["10K_norm"]

# Location Loop
for location in locationList:

  # Matrix List
  matrixList = ["3B9_5-.txt", "3B9_5plus.txt", "69_127-.txt", "69_127plus.txt"]

  # Matrix Loop
  for matrix in matrixList:

    # Region files
    regionList = ["STAG1_only.bed", "STAG2_only.bed"]

    # Loop on region files
    for region in regionList:

      # Parameters
      if("10K_" in location):
        res = "10000"
        minDist = "100000"
      elif("25K_" in location):
        res = "25000"
        minDist = "250000"
      elif("50K_" in location):
        res = "50000"
        minDist = "500000"
      matrixName = il+location+"/"+matrix
      regionName = rl+region
      outName = location+"/"+matrix+":"+region.split("_")[0]
      print outName

      # Input
      loopBins = "10"
      resolution = res
      medThresh = "25"
      highThresh = "70"
      binsMiddle = "0"
      acceptedErrors = "0"
      minimumDistance = minDist
      thresholdFileName = matrixName[:-4]+"_percentiles.txt"
      regionFileName = regionName
      matrixFileName = matrixName
      outputRegionsFileName = olm+outName+"_regions.txt"
      outputMatrixFileName = olm+outName+".txt"
      outputGraphFileName = olg+outName+".pdf"

      # Creating folders
      f1 = "/".join(outputRegionsFileName.split("/")[:-1])+"/"
      os.system("mkdir -p "+f1)
      f2 = "/".join(outputGraphFileName.split("/")[:-1])+"/"
      os.system("mkdir -p "+f2)

      # Creating matrix
      command = "python loopOgram.py "+" ".join([loopBins, resolution, medThresh, highThresh, binsMiddle, acceptedErrors, minimumDistance, thresholdFileName, regionFileName, matrixFileName, outputRegionsFileName, outputMatrixFileName])
      os.system(command)

      # Creating heatmap
      #command = "R CMD BATCH '--args '"+outputMatrixFileName+"' '"+outputRegionsFileName+"' '"+outputGraphFileName+" loopOgramHeatmap.R loopOgramHeatmap.Rout"
      #os.system(command)


