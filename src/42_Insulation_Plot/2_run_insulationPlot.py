
# Import
import os
import sys
from glob import glob

# Stag List
il = "/home/egg/Projects/Wendt_Stag/Results/11_Insulation_Plot/1_Insulation_Table/25K_norm/"
ol = "/home/egg/Projects/Wendt_Stag/Results/11_Insulation_Plot/2_Insulation_Plots/white_orange_red/25K_norm/"
tableList = glob(il + "*/*.txt")

# Stag Loop
for tableFileName in tableList:

  # Parameters
  if("minus" in tableFileName): foldChange = "1"
  else: foldChange = "0"
  outLoc = ol + tableFileName.split("/")[-2] + "/"
  os.system("mkdir -p "+outLoc)
  fileName = tableFileName.split("/")[-1].split(".")[0] + ".pdf"

  # Input
  foldChange = foldChange
  regionName = tableFileName.split("/")[-2]
  signalName = tableFileName.split("/")[-1].split(".")[0]
  inputMatrixFileName = tableFileName
  outputFileName = outLoc + fileName

  # Write
  command = "Rscript 2_insulationPlot.R "+" ".join([foldChange, regionName, signalName, inputMatrixFileName, outputFileName])
  os.system(command)


