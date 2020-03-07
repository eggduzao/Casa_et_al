
# Import
import os
import sys
from glob import glob

# Stag List
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/40_Fig4G_Insulation_RaoNeg/1_Insulation_Table/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/40_Fig4G_Insulation_RaoNeg/2_Insulation_Plots/"
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


