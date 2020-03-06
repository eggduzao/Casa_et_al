
# Import
import os
import sys

# Condition List
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/23_AB_Compartment_Losada/2_Compartments/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/23_AB_Compartment_Losada/3_AB_Plot/"
condList = [["SISA1", "LCONT"], ["SISA2", "LCONT"]]

# Condition Loop
for cond in condList:

  # Input
  abTreatFileName = il + cond[0] + ".bed"
  abContrFileName = il + cond[1] + ".bed"
  command = "cat " + il + cond[0] + "/*.bed > " + il + cond[0] + ".bed"
  os.system(command)
  command = "cat " + il + cond[1] + "/*.bed > " + il + cond[1] + ".bed"
  os.system(command)
  outputFileName = ol + cond[0] + "_AB.txt"
  outputBedFilePrefix = ol + cond[0] + "_AB"

  # Execution
  command = "python 3_createPlotTable.py "+" ".join([abTreatFileName, abContrFileName, outputFileName, outputBedFilePrefix])
  os.system(command)


