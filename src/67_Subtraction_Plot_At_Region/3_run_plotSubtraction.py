
# Import
import os
import sys
from glob import glob

# Matrix List
itisstd = "N"
regionFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/0_Regions/Selected_Regions2.txt"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/input/"
ml1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/1_Bin_Matrices/selected_regions_2/"
ml2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/2_Subtraction_Matrices/selected_regions_2/"
olm1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/3_Matrix_Visualization/selected_regions_2/1_sub_matrix/"
olp1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/3_Matrix_Visualization/selected_regions_2/1_sub_plot_png/"
olm2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/3_Matrix_Visualization/selected_regions_2/2_org_matrix/"
olp2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/38_Subtraction_Plot_At_Region/3_Matrix_Visualization/selected_regions_2/2_org_plot_png/"

###################################################################################################
# SUBTRACTED
###################################################################################################

# Opening files
inputFile1Name = fl + "31_rps.txt"
inputFile2Name = fl + "32_rps.txt"
inputFile1 = open(inputFile1Name, "w")
inputFile2 = open(inputFile2Name, "w")

# Region Loop
regionFile = open(regionFileName, "rU")
for line in regionFile:

  ll = line.strip().split("\t")

  # Hic Name List
  hicList = ["STAG1", "STAG2"]

  # Hic Name Loop
  for hicname in hicList:

    # Name
    fullnamebin = "_".join([hicname, ll[0], ll[1], ll[2], "bin"])
    fullnameraw = "_".join([hicname, ll[0], ll[1], ll[2], "raw"])

    # Condition List
    condList = ["25K_norm"]

    # Condition Loop
    for cond in condList:
 
      # Parameters
      if("25K" in cond): res = 25000
      elif("100K" in cond): res = 100000

      # Percentile List
      percList = ["80", "90"]

      # Percentile Loop
      for perc in percList:

        # Input
        chrom = ll[0]
        pos1 = ll[1]
        pos2 = ll[2]
        resolution = str(res)
        minV = "-1"
        maxV = "1"
        matrixFileName = ml2 + perc + "/" + cond + "/" + fullnamebin + ".txt"
        outputMatrixFileName = olm1 + perc + "/" + cond + "/" + fullnamebin + ".txt"
        outputPlotFileNameAll = olp1 + perc + "/" + cond + "/" + fullnamebin + "_ALL.png"
        outputPlotFileNameBlue = olp1 + perc + "/" + cond + "/" + fullnamebin + "_DECREASE.png"
        outputPlotFileNameRed = olp1 + perc + "/" + cond + "/" + fullnamebin + "_INCREASE.png"
        os.system("mkdir -p " + olp1 + perc + "/" + cond + "/")

        # Creating matrix
        inputFile1.write(" ".join([chrom, pos1, pos2, resolution, matrixFileName, outputMatrixFileName])+"\n")

        # Creating plot
        inputFile2.write(" ".join([itisstd, minV, maxV, outputMatrixFileName, outputPlotFileNameAll, outputPlotFileNameBlue, outputPlotFileNameRed])+"\n")

        # Range List
        minVList = [-1, -5, -10, -20, -30, -50, -100]
        maxVList = [1, 5, 10, 20, 30, 50, 100]

        # Range Loop
        for i in range(0,len(minVList)):

          newfullnameraw = "_".join([fullnameraw, str(minVList[i]), str(maxVList[i])])

          # Input
          chrom = ll[0]
          pos1 = ll[1]
          pos2 = ll[2]
          resolution = str(res)
          minV = str(minVList[i])
          maxV = str(maxVList[i])
          matrixFileName = ml2 + perc + "/" + cond + "/" + fullnameraw + ".txt"
          outputMatrixFileName = olm1 + perc + "/" + cond + "/" + newfullnameraw + ".txt"
          outputPlotFileNameAll = olp1 + perc + "/" + cond + "/" + newfullnameraw + "_ALL.png"
          outputPlotFileNameBlue = olp1 + perc + "/" + cond + "/" + newfullnameraw + "_DECREASE.png"
          outputPlotFileNameRed = olp1 + perc + "/" + cond + "/" + newfullnameraw + "_INCREASE.png"
          os.system("mkdir -p " + olp1 + perc + "/" + cond + "/")

          # Creating matrix
          inputFile1.write(" ".join([chrom, pos1, pos2, resolution, matrixFileName, outputMatrixFileName])+"\n")

          # Creating plot
          inputFile2.write(" ".join([itisstd, minV, maxV, outputMatrixFileName, outputPlotFileNameAll, outputPlotFileNameBlue, outputPlotFileNameRed])+"\n")

###################################################################################################
# NOT SUBTRACTED
###################################################################################################

# Region Loop
regionFile = open(regionFileName, "rU")
for line in regionFile:

  ll = line.strip().split("\t")

  # Hic Name List
  hicList = ["STAG1_P", "STAG1_M", "STAG2_P", "STAG2_M"]

  # Hic Name Loop
  for hicName in hicList:

    hicname = "_".join([hicName, ll[0], ll[1], ll[2], "raw"])

    # Condition List
    condList = ["25K_norm"]

    # Condition Loop
    for cond in condList:
 
      # Parameters
      if("25K" in cond): res = 25000
      elif("100K" in cond): res = 100000

      # Percentile List
      percList = ["80", "90"]

      # Percentile Loop
      for perc in percList:

        # Range List
        minVList = [-1, -5, -10, -20, -30, -50, -100]
        maxVList = [1, 5, 10, 20, 30, 50, 100]

        # Range Loop
        for i in range(0,len(minVList)):

          newfullnameraw = "_".join([hicname, str(minVList[i]), str(maxVList[i])])

          # Input
          chrom = ll[0]
          pos1 = ll[1]
          pos2 = ll[2]
          resolution = str(res)
          minV = str(minVList[i])
          maxV = str(maxVList[i])
          matrixFileName = ml1 + perc + "/" + cond + "/" + hicname + ".txt"
          outputMatrixFileName = olm2 + perc + "/" + cond + "/" + newfullnameraw + ".txt"
          outputPlotFileNameAll = olp2 + perc + "/" + cond + "/" + newfullnameraw + "_ALL.png"
          outputPlotFileNameBlue = olp2 + perc + "/" + cond + "/" + newfullnameraw + "_DECREASE.png"
          outputPlotFileNameRed = olp2 + perc + "/" + cond + "/" + newfullnameraw + "_INCREASE.png"
          os.system("mkdir -p " + olp2 + perc + "/" + cond + "/")

          # Creating matrix
          inputFile1.write(" ".join([chrom, pos1, pos2, resolution, matrixFileName, outputMatrixFileName])+"\n")

          # Creating plot
          inputFile2.write(" ".join([itisstd, minV, maxV, outputMatrixFileName, outputPlotFileNameAll, outputPlotFileNameBlue, outputPlotFileNameRed])+"\n")

inputFile1.close()
inputFile2.close()


