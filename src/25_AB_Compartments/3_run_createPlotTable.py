
# Import
import os
import sys

# Condition List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/25_AB_Compartments/2_compartments/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/25_AB_Compartments/3_conservation_plot/"
condList = ["STAG1", "STAG2"]

# Condition Loop
for cond in condList:

  # Input
  abTreatFileName = il+cond+"+.bed"
  abContrFileName = il+cond+"-.bed"
  outputFileName = ol+cond+".txt"
  outputBedFilePrefix = ol+cond

  # Execution
  command = "python 3_createPlotTable.py "+" ".join([abTreatFileName, abContrFileName, outputFileName, outputBedFilePrefix])
  os.system(command)

# Input
abTreatFileName = il+"STAG1-.bed"
abContrFileName = il+"STAG2-.bed"
outputFileName = ol+"WT.txt"
outputBedFilePrefix = ol+"WT"

# Execution
command = "python 3_createPlotTable.py "+" ".join([abTreatFileName, abContrFileName, outputFileName, outputBedFilePrefix])
os.system(command)


