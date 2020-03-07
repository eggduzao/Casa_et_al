
# Import
import os
import sys
from glob import glob

# Stag List
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/40_Fig4G_Insulation_RaoNeg/2_Insulation_Plots/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/40_Fig4G_Insulation_RaoNeg/3_Insulation_Plots_2/"
tableList = glob(il + "*/*.tsv")

# Stag Loop
for tableFileName in tableList:

  # Parameters
  if("minus" in tableFileName): fcc = "1"
  else: fcc = "0"
  outLoc = ol + tableFileName.split("/")[-2] + "/"
  os.system("mkdir -p "+outLoc)
  fileName = tableFileName.split("/")[-1].split(".")[0] + ".pdf"

  # Input
  foldChange = fcc
  inputMatrixFileName = tableFileName
  outputFileName = outLoc + fileName

  # Write
  command = "Rscript 3_insulationPlot2.R "+" ".join([foldChange, inputMatrixFileName, outputFileName])
  os.system(command)


