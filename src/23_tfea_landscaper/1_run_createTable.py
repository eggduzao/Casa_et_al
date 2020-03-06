
# Import
import os
import sys

# Tfea List
tl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/23_tfea_landscaper/0_raw_tfea/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/23_tfea_landscaper/1_table_tfea/"
tfeaList = [["R_STAG1_RM_STAG1_intersection_Betas.csv", "R_STAG1_RM_STAG1_minusOnly_Betas.csv"], ["R_STAG1_RM_STAG1_intersection_Betas.csv", "R_STAG1_RM_STAG1_plusOnly_Betas.csv"], ["R_STAG1_RM_STAG1_minusOnly_Betas.csv", "R_STAG1_RM_STAG1_plusOnly_Betas.csv"], ["STAG1only_mrg_replicates_Betas.csv", "STAG2only_mtg_replicates_Betas.csv"]]
labelList = [["intersection", "minus"], ["intersection", "plus"], ["minus", "plus"], ["STAG1", "STAG2"]]

# Tfea Loop
for i in range(0,len(tfeaList)):
  
  # Parameters
  name = "_".join(labelList[i])

  # Input
  inputFileName1 = tl+tfeaList[i][0]
  inputFileName2 = tl+tfeaList[i][1]
  outputFileName = ol+name+".txt"

  # Execution
  command = "python 1_createTable.py "+" ".join([inputFileName1, inputFileName2, outputFileName])
  os.system(command)


