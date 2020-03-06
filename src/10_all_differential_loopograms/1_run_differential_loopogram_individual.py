
# Import
import os
import sys

# Location List
ilm = "/home/egg/Projects/Wendt_Stag/Data/stag_matrix_files/"
ilr = "/home/egg/Projects/Wendt_Stag/Previous_Results/3_Loopograms/matrix_new/"
olm = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/individual_test/matrix/"
olg = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/individual_test/graphs/"
locationList = ["10K_none", "10K_norm", "25K_none", "25K_norm", "50K_none", "50K_norm"]
locationList = ["10K_norm"]

# Location Loop
for location in locationList:

  # Matrix List
  matrixList = ["3B9_5-.txt", "3B9_5plus.txt", "69_127-.txt", "69_127plus.txt"]

  # Matrix Loop
  for matrix in matrixList:

    # Region files
    regionList = ["3B9_5-.txt:STAG1_regions.txt", "3B9_5-.txt:STAG2_regions.txt", "3B9_5plus.txt:STAG1_regions.txt", 
                  "3B9_5plus.txt:STAG2_regions.txt", "69_127-.txt:STAG1_regions.txt", "69_127-.txt:STAG2_regions.txt",
                  "69_127plus.txt:STAG1_regions.txt", "69_127plus.txt:STAG2_regions.txt"]

    # Loop on region files
    for region in regionList:

      # Parameters
      if(matrix == "3B9_5-.txt"): mLabel = "M_STAG2minus"
      elif(matrix == "3B9_5plus.txt"): mLabel = "M_STAG2plus"
      elif(matrix == "69_127-.txt"): mLabel = "M_STAG1minus"
      elif(matrix == "69_127plus.txt"): mLabel = "M_STAG1plus"

      if("3B9_5-.txt" in region): rmLabel = "RM_STAG2minus"
      elif("3B9_5plus.txt" in region): rmLabel = "RM_STAG2plus"
      elif("69_127-.txt" in region): rmLabel = "RM_STAG1minus"
      elif("69_127plus.txt" in region): rmLabel = "RM_STAG1plus"

      if("STAG1" in region): rLabel = "R_STAG1"
      elif("STAG2" in region): rLabel = "R_STAG2"

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
    

