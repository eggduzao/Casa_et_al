
# Import
import os
import sys

# Location List
ilm = "/home/egg/Projects/Wendt_Stag/Data/stag_matrix_files/"
ilr = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/input/"
olm = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/combined_test/matrix/"
olg = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/combined_test/graphs/"
locationList = ["10K_none", "10K_norm", "25K_none", "25K_norm", "50K_none", "50K_norm"]
locationList = ["10K_norm"]

# Location Loop
for location in locationList:

  # Matrix List
  matrixList = ["3B9_5-.txt", "3B9_5plus.txt", "69_127-.txt", "69_127plus.txt"]

  # Matrix Loop
  for matrix in matrixList:

    # Region files
    regionList = ["R_STAG1_RM_STAG1_intersection.bed", "R_STAG1_RM_STAG1_minusOnly.bed", "R_STAG1_RM_STAG1_plusOnly.bed"]

    # Loop on region files
    for region in regionList:

      # Parameters
      if(matrix == "3B9_5-.txt"): mLabel = "M_STAG2minus"
      elif(matrix == "3B9_5plus.txt"): mLabel = "M_STAG2plus"
      elif(matrix == "69_127-.txt"): mLabel = "M_STAG1minus"
      elif(matrix == "69_127plus.txt"): mLabel = "M_STAG1plus"

      rr = region.split(".")[0].split("_")
      rLabel = "_".join(rr[:2])
      rmLabel = "_".join(rr[2:])

      name = "_".join([rLabel,rmLabel,mLabel])

      if("10K_" in location): res = "10000"
      elif("25K_" in location): res = "25000"
      elif("50K_" in location): res = "50000"

      # Input
      loopBins = "10"
      resolution = res
      regionFileName = ilr+location+"/"+region
      matrixFileName = ilm+location+"/"+matrix
      outputMatrixFileName = olm+location+"/"+name+".txt"
      outputGraphFileName = olg+location+"/"+name+".pdf"
      os.system("mkdir -p "+olm+location)
      os.system("mkdir -p "+olg+location)

      # Creating matrix
      #command = "python 1_differential_loopogram.py "+" ".join([loopBins, resolution, regionFileName, matrixFileName, outputMatrixFileName])
      #os.system(command)

      # Creating heatmap
      command = "R CMD BATCH '--args '"+outputMatrixFileName+"' '"+regionFileName+"' '"+outputGraphFileName+" 1_loopOgramHeatmap.R 1_loopOgramHeatmap.Rout"
      os.system(command)
    

