
# Import
import os
import sys

# Gene List
il1 = "/home/egg/Projects/Wendt_Stag/Results/11_Insulation_Plot/0_Input_Regions/"
il2 = "/home/egg/Projects/Wendt_Stag/Data/stag_bed_files/"
il3 = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/input/10K_norm/"
ml = "/home/egg/Projects/Wendt_Stag/Data/stag_matrix_files/10K_norm/"
ol = "/home/egg/Projects/Wendt_Stag/Results/11_Insulation_Plot/1_Insulation_Table/10K_norm/"
geneList = ["STAG1_downregulated_TSS", "STAG1_upregulated_TSS", "STAG2_downregulated_TSS", "STAG2_upregulated_TSS",
            "STAG1_downregulated_TTS", "STAG1_upregulated_TTS", "STAG2_downregulated_TTS", "STAG2_upregulated_TTS",
            "STAG1_only", "STAG2_only", "shared",
            "R_STAG1_RM_STAG1_intersection", "R_STAG1_RM_STAG1_minusOnly", "R_STAG1_RM_STAG1_plusOnly", "R_STAG2_RM_STAG2_intersection", "R_STAG2_RM_STAG2_minusOnly", "R_STAG2_RM_STAG2_plusOnly"]

# Gene Loop
for geneName in geneList:

  # Matrix List
  if("STAG1" in geneName):
    matrixListList = [["69_127-"], ["69_127plus"], ["69_127plus", "69_127-"]]
    labelList = ["STAG1_WT", "STAG1_DEG", "STAG1_DEG_minus_STAG1_WT"]
  elif("STAG2" in geneName):
    matrixListList = [["3B9_5-"], ["3B9_5plus"], ["3B9_5plus", "3B9_5-"]]
    labelList = ["STAG2_WT", "STAG2_DEG", "STAG2_DEG_minus_STAG2_WT"]
  else:
    matrixListList = [["69_127-"], ["69_127plus"], ["69_127plus", "69_127-"], ["3B9_5-"], ["3B9_5plus"], ["3B9_5plus", "3B9_5-"]]
    labelList = ["STAG1_WT", "STAG1_DEG", "STAG1_DEG_minus_STAG1_WT", "STAG2_WT", "STAG2_DEG", "STAG2_DEG_minus_STAG2_WT"]

  il = il1
  if("only" in geneName or "shared" in geneName): il = il2
  elif("RM" in geneName): il = il3

  # Stag Loop
  for i in range(0,len(matrixListList)):

    # Parameters
    matrixList = matrixListList[i]
    label = labelList[i]

    # Input
    halfBin = "20"
    resolution = "10000"
    regionFileName = il + geneName + ".bed"
    matrix1FileName = ml + matrixList[0] + ".txt"
    if(len(matrixList) > 1): matrix2FileName = ml + matrixList[1] + ".txt"
    else: matrix2FileName = "NA"
    outputTableFileName = ol + geneName + "/" + label + ".txt"

    # Write
    command = "python 1_insulationTable.py "+" ".join([halfBin, resolution, regionFileName, matrix1FileName, matrix2FileName, outputTableFileName])
    os.system(command)


