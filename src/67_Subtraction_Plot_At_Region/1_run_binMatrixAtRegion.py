
# Import
import os
import sys

# Region List
regionFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/0_Regions/Selected_Regions2.txt"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/"
tl = "/scratch/egadegu/CBF/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/1_Bin_Matrices/selected_regions_2/"

# Open file
inputFileName = fl + "1_cbf.txt"
inFile = open(inputFileName, "w")

# Region Loop
regionFile = open(regionFileName, "rU")
for line in regionFile:

  ll = line.strip().split("\t")

  # Hic List
  hicNameList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

  # Hic Loop
  for hicName in hicNameList:

    # Condition List
    condList = ["25K_norm"]

    # Condition Loop
    for cond in condList:

      # Parameters
      if("25K" in cond): res = 25000
      elif("100K" in cond): res = 100000
      if(hicName == "69_127plus"): hicnewname = "STAG1_P"
      elif(hicName == "69_127-"): hicnewname = "STAG1_M"
      elif(hicName == "3B9_5plus"): hicnewname = "STAG2_P"
      elif(hicName == "3B9_5-"): hicnewname = "STAG2_M"
      outName = "_".join([hicnewname, ll[0], ll[1], ll[2]])

      # Percentile List
      percList = ["80", "90"]

      # Percentile Loop
      for perc in percList:

        # Input
        chromosome = ll[0]
        region1 = ll[1]
        region2 = ll[2]
        resolution = str(res)
        percentileThreshold = perc
        percentileFileName = il + cond + "/" + hicName + "_percentiles.txt"
        inputMatrixFileName = il + cond + "/" + hicName + ".txt"
        tempLocation = tl + perc + "/" + cond + "/" + outName + "/"
        outBinPrefix = ol + perc + "/" + cond + "/" + outName + "_bin"
        outRawPrefix = ol + perc + "/" + cond + "/" + outName + "_raw"

        # Write
        inFile.write(" ".join([chromosome, region1, region2, resolution, percentileThreshold, percentileFileName, inputMatrixFileName, tempLocation, outBinPrefix, outRawPrefix])+"\n")

inFile.close()


