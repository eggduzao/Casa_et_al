
###################################################################################################
# Input + Execution
###################################################################################################

# Import
import os
import sys

#################################################
# BIN / STD
#################################################

# Matrix List
chromSizesFile = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/input/"
ml = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/2_Subtraction_Matrices/"
olm = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/4_Matrix_Visualization/matrix/"
olp = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/4_Matrix_Visualization/plot_png/"
matrixList = ["3B9_5plus_bin_M_3B9_5-_bin", "3B9_5plus_std_M_3B9_5-_std", "69_127plus_bin_M_69_127-_bin", "69_127plus_std_M_69_127-_std"] 

# Opening files
inputFile1Name = fl + "31_fsp.txt"
inputFile2Name = fl + "32_fsp.txt"
inputFile1 = open(inputFile1Name, "w")
inputFile2 = open(inputFile2Name, "w")

# Matrix Loop
for matrix in matrixList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23) + ["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Condition List
    condList = ["25K_none", "25K_norm", "50K_none", "50K_norm", "250K_none", "250K_norm"]
    #condList = ["100K_none", "100K_norm"]

    # Condition Loop
    for cond in condList:

      # Parameters
      if("_bin_" in matrix):
        itisstd = "N"
        minV = "-1"
        maxV = "1"
      else:
        itisstd = "Y"
        minV = "-5"
        maxV = "5"

      # Parameters
      if("10K" in cond): res = 10000
      elif("25K" in cond): res = 25000
      elif("50K" in cond): res = 50000
      elif("100K" in cond): res = 100000
      elif("250K" in cond): res = 250000
      elif("500K" in cond): res = 500000
      elif("1000K" in cond): res = 1000000
 
      # Parameters
      if(matrix == "3B9_5plus_bin_M_3B9_5-_bin"): outName = "STAG2_bin"
      elif(matrix == "3B9_5plus_std_M_3B9_5-_std"): outName = "STAG2_std"
      elif(matrix == "69_127plus_bin_M_69_127-_bin"): outName = "STAG1_bin"
      elif(matrix == "69_127plus_std_M_69_127-_std"): outName = "STAG1_std"

      # Percentile List
      percList = ["70", "80", "90", "95", "98"]

      # Percentile Loop
      for perc in percList:

        # Input
        chrom = chrom
        resolution = str(res)
        chromSizesFileName = chromSizesFile
        matrixFileName = ml + perc + "/" + cond + "/" + matrix + ".txt"
        outputMatrixFileName = olm + perc + "/" + cond + "/" + matrix + "_" + chrom + ".txt"
        outputPlotFileNameAll = olp + perc + "/" + cond + "/" + outName + "_" + chrom + "_ALL.png"
        outputPlotFileNameBlue = olp + perc + "/" + cond + "/" + outName + "_" + chrom + "_DECREASE.png"
        outputPlotFileNameRed = olp + perc + "/" + cond + "/" + outName + "_" + chrom + "_INCREASE.png"
        os.system("mkdir -p " + olp + perc + "/" + cond + "/")

        # Creating matrix
        inputFile1.write(" ".join([chrom, resolution, chromSizesFileName, matrixFileName, outputMatrixFileName])+"\n")

        # Creating plot
        inputFile2.write(" ".join([itisstd, minV, maxV, outputMatrixFileName, outputPlotFileNameAll, outputPlotFileNameBlue, outputPlotFileNameRed])+"\n")

#################################################
# RAW
#################################################

# Matrix List
matrixList = ["3B9_5plus_M_3B9_5-", "69_127plus_M_69_127-"] 

# Matrix Loop
for matrix in matrixList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23) + ["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Condition List
    condList = ["25K_none", "25K_norm", "50K_none", "50K_norm", "250K_none", "250K_norm"]
    #condList = ["100K_none", "100K_norm"]

    # Condition Loop
    for cond in condList:

      # Parameters
      itisstd = "Y"
      minV = "-30"
      maxV = "30"

      # Parameters
      if("10K" in cond): res = 10000
      elif("25K" in cond): res = 25000
      elif("50K" in cond): res = 50000
      elif("100K" in cond): res = 100000
      elif("250K" in cond): res = 250000
      elif("500K" in cond): res = 500000
      elif("1000K" in cond): res = 1000000
 
      # Parameters
      if(matrix == "3B9_5plus_M_3B9_5-"): outName = "STAG2_raw"
      elif(matrix == "69_127plus_M_69_127-"): outName = "STAG1_raw"

      # Percentile List
      percList = ["raw"]

      # Percentile Loop
      for perc in percList:

        # Input
        chrom = chrom
        resolution = str(res)
        chromSizesFileName = chromSizesFile
        matrixFileName = ml + perc + "/" + cond + "/" + matrix + ".txt"
        outputMatrixFileName = olm + perc + "/" + cond + "/" + matrix + "_" + chrom + ".txt"
        outputPlotFileNameAll = olp + perc + "/" + cond + "/" + outName + "_" + chrom + "_ALL.png"
        outputPlotFileNameBlue = olp + perc + "/" + cond + "/" + outName + "_" + chrom + "_DECREASE.png"
        outputPlotFileNameRed = olp + perc + "/" + cond + "/" + outName + "_" + chrom + "_INCREASE.png"
        os.system("mkdir -p " + olp + perc + "/" + cond + "/")

        # Creating matrix
        inputFile1.write(" ".join([chrom, resolution, chromSizesFileName, matrixFileName, outputMatrixFileName])+"\n")

        # Creating plot
        inputFile2.write(" ".join([itisstd, minV, maxV, outputMatrixFileName, outputPlotFileNameAll, outputPlotFileNameBlue, outputPlotFileNameRed])+"\n")

inputFile1.close()
inputFile2.close()


