
# Import
import os
import sys

# Condition List
counter = 1
fl = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/5_Loops/input/"
rl = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/5_Loops/3_contact_files_filtered/"
ml = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_matrix_files/25K_norm/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/5_Loops/5_ChoosingBestLoops/LoopOGrams/"
condNameList = [["69_127-", "3B9_5-"]]
condLabelList = ["STAG12_FC"]

# Condition Loop
for i in range(0,len(condNameList)):

  # Parameters
  name1 = condNameList[i][0]
  name2 = condNameList[i][1]
  label = condLabelList[i]

  # Region List
  regionList = ["STAG1_minus_auxin_50Kb_5Mb", "STAG1_minus_auxin_50Kb_10Mb", "STAG1_minus_auxin_150Kb_5Mb" , "STAG1_minus_auxin_150Kb_10Mb",
                "STAG1_plus_auxin_50Kb_5Mb", "STAG1_plus_auxin_50Kb_10Mb", "STAG1_plus_auxin_150Kb_5Mb" , "STAG1_plus_auxin_150Kb_10Mb",
                "STAG2_minus_auxin_50Kb_5Mb", "STAG2_minus_auxin_50Kb_10Mb", "STAG2_minus_auxin_150Kb_5Mb" , "STAG2_minus_auxin_150Kb_10Mb",
                "STAG2_plus_auxin_50Kb_5Mb", "STAG2_plus_auxin_50Kb_10Mb", "STAG2_plus_auxin_150Kb_5Mb" , "STAG2_plus_auxin_150Kb_10Mb"]
  
  # Region Loop
  for region in regionList:

    # Parameters
    outLoc = ol + "signal_" + label + "/"
    os.system("mkdir -p "+outLoc)

    # Input
    loopBins = "10"
    resolution = "25000"
    regionsFileName = rl + region + ".txt"
    hicFile1Name = ml + name1 + ".txt"
    hicFile2Name = ml + name2 + ".txt"
    outputFileName = outLoc + region + ".txt"

    # Write
    inFile = open(fl + str(counter) + "_cbl.txt", "w")
    inFile.write("\n".join([loopBins, resolution, regionsFileName, hicFile1Name, hicFile2Name, outputFileName]))
    inFile.close()
    counter += 1


