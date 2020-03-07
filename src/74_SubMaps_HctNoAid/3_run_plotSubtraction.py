
# Import
import os
import sys

# Matrix List
itisstd = "N"
res = "25000"
regionFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/0_Regions/Selected_Regions.txt"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/"
ml1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/2_Subtraction_Matrices/"
ml2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/1_Bin_Matrices/"
olm1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/3_Matrix_Visualization/1_sub_matrix/"
olp1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/3_Matrix_Visualization/1_sub_plot_png/"
olm2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/3_Matrix_Visualization/2_org_matrix/"
olp2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/45_FigS8B_SubMaps_RaoNeg/3_Matrix_Visualization/2_org_plot_png/"

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
  hicList = ["STAG1+_v_Untreated", "STAG2+_v_Untreated"]

  # Hic Name Loop
  for hicname in hicList:

    # Parameters
    fullnameraw = "_".join([hicname, ll[0], ll[1], ll[2], "raw"])
    newfullnameraw = "_".join([fullnameraw, "-10", "10"])

    # Input
    chrom = ll[0]
    pos1 = ll[1]
    pos2 = ll[2]
    resolution = res
    minV = "-10"
    maxV = "10"
    matrixFileName = ml1 + fullnameraw + ".txt"
    outputMatrixFileName = olm1 + newfullnameraw + ".txt"
    outputPlotFileNameAll = olp1 + newfullnameraw + "_ALL.png"
    outputPlotFileNameBlue = olp1 + newfullnameraw + "_DECREASE.png"
    outputPlotFileNameRed = olp1 + newfullnameraw + "_INCREASE.png"

    # Creating matrix
    inputFile1.write(" ".join([chrom, pos1, pos2, resolution, matrixFileName, outputMatrixFileName])+"\n")

    # Creating plot
    inputFile2.write(" ".join([itisstd, minV, maxV, outputMatrixFileName, outputPlotFileNameAll, outputPlotFileNameBlue, outputPlotFileNameRed])+"\n")

regionFile.close()

###################################################################################################
# NOT SUBTRACTED
###################################################################################################

# Region Loop
regionFile = open(regionFileName, "rU")
for line in regionFile:

  ll = line.strip().split("\t")

  # Hic Name List
  hicList = ["STAG1+", "STAG2+", "Untreated"]

  # Hic Name Loop
  for hicName in hicList:

    hicname = "_".join([hicName, ll[0], ll[1], ll[2], "raw"])
    newfullnameraw = "_".join([hicname, "-10", "10"])

    # Input
    chrom = ll[0]
    pos1 = ll[1]
    pos2 = ll[2]
    resolution = res
    minV = "-10"
    maxV = "10"
    matrixFileName = ml2 + hicname + ".txt"
    outputMatrixFileName = olm2 + newfullnameraw + ".txt"
    outputPlotFileNameAll = olp2 + newfullnameraw + "_ALL.png"
    outputPlotFileNameBlue = olp2 + newfullnameraw + "_DECREASE.png"
    outputPlotFileNameRed = olp2 + newfullnameraw + "_INCREASE.png"

    # Creating matrix
    inputFile1.write(" ".join([chrom, pos1, pos2, resolution, matrixFileName, outputMatrixFileName])+"\n")

    # Creating plot
    inputFile2.write(" ".join([itisstd, minV, maxV, outputMatrixFileName, outputPlotFileNameAll, outputPlotFileNameBlue, outputPlotFileNameRed])+"\n")

inputFile1.close()
inputFile2.close()
regionFile.close()


