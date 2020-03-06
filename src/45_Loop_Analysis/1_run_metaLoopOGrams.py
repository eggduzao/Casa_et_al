
# Import
import os
import sys

# Condition List
counter = 1
fl = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/8_Meta_Plots/input/"
rl = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/8_Meta_Plots/0_Regions/"
ml = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_matrix_files/25K_norm/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/8_Meta_Plots/1_Meta_LoopOGram/"
condNameList = [["69_127-", "69_127plus"], ["3B9_5-", "3B9_5plus"]]
condLabelList = ["STAG1_FC", "STAG2_FC"]

# Condition Loop
for i in range(0,len(condNameList)):

  # Parameters
  name1 = condNameList[i][0]
  name2 = condNameList[i][1]
  label = condLabelList[i]

  # Region List
  regionList1 = ["STAG1_up", "STAG2_up"]
  regionList2 = ["STAG1_down", "STAG2_down"]
  
  # Region Loop
  for j in range(0,len(regionList1)):

    # Parameters
    region1 = regionList1[j]
    region2 = regionList2[j]
    outLoc = ol + "region_" + region1.split("_")[0] + "_up_down" + "/"
    os.system("mkdir -p "+outLoc)

    # Input
    loopBins = "10"
    resolution = "25000"
    regionsFile1Name = rl + region1 + ".bed"
    regionsFile2Name = rl + region2 + ".bed"
    hicFile1Name = ml + name1 + ".txt"
    hicFile2Name = ml + name2 + ".txt"
    outputFileName = outLoc + label + ".txt"

    # Write
    inFile = open(fl + str(counter) + "_mlg.txt", "w")
    inFile.write("\n".join([loopBins, resolution, regionsFile1Name, regionsFile2Name, hicFile1Name, hicFile2Name, outputFileName]))
    inFile.close()
    counter += 1


