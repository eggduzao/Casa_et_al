
# Import
import os
import sys

# Bam List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/2_Merged_Bam_Files/"
tl = "/scratch/egadegu/MAC/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/4_Peaks/"
treatList = ["Smc3_Treated", "Smc3_Untreated"]
controlList = ["Smc3_Input", "Smc3_Input"]

# Opening file
inFileName = fl + "4_mac.txt"
inFile = open(inFileName, "w")

# Bam Loop
for i in range(0,len(treatList)):

  # Parameters
  treatName = treatList[i]
  controlName = controlList[i]

  # Analysis Type List
  anTypeList = ["Ed", "Macarena"]

  # Analysis Type Loop
  for anType in anTypeList:

    # Pvalue List
    pvalueList = ["0.05", "0.01", "0.005", "0.001", "0.0005", "0.0001"]

    # Pvalue Loop
    for pvalue in pvalueList:

      # Parameters
      pv = pvalue.split(".")[-1]
      outName = "_".join([treatName, anType, pv])

      # Input
      name = outName
      pValue = pvalue
      commandType = anType
      treatmentFileName = il + treatName + ".bam"
      controlFileName = il + controlName + ".bam"
      tempLocation = tl + outName + "/"
      outputLocation = ol

      # Creating files
      inFile.write(" ".join([name, pValue, commandType, treatmentFileName, controlFileName, tempLocation, outputLocation])+"\n")

# Close
inFile.close()


