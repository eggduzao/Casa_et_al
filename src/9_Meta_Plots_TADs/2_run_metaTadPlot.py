# Import
import os
import sys

# Stag List
counter = 1
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/8_Meta_Plots/input/"
il = "/projects/ag-papan/eduardo/Wendt_Stag/Results/8_Meta_Plots/4_Meta_Tad_Tables/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/8_Meta_Plots/5_Meta_Tad_Plots/"
stagList = ["STAG1", "STAG2", "WT"]

# Stag Loop
for stag in stagList:

  # Diff Type List
  diffTypeList = ["0_equal", "1_split", "2_merge", "3_shift_upstream", "4_shift_downstream", "5_new"]

  # Diff Type Loop
  for diff in diffTypeList:

    # Signal List
    signalLabelList = [["STAG1+","STAG1-"], ["STAG2+","STAG2-"], ["STAG1-","STAG2-"], ["STAG1+","STAG2+"]]

    # Signal Loop
    for signalLabel in signalLabelList:

      # Parameters
      signalLabel = "_vs_".join(signalLabel)
      inLoc = il + "/".join([stag, diff]) + "/"
      outLoc = ol + "/".join([stag, diff]) + "/"
      command = "mkdir -p "+outLoc
      os.system(command)

      # Input
      inputTableFileName = inLoc + signalLabel + ".txt"
      outputFileName = inLoc + signalLabel + ".pdf"

      # Write
      inFile = open(fl + str(counter) + "_mtp.txt", "w")
      inFile.write("\n".join([inputTableFileName, outputFileName]))
      inFile.close()
      counter += 1


