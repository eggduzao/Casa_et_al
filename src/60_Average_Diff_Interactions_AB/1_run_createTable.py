
# Import
import os
import sys

# List
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/31_Average_Diff_Interactions_AB/input/"
al = "/usr/users/egadegu/Projects/Wendt_Stag/Results/3_AB_Compartments/2_compartments/"
tl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/4_TADs/1_TADs_GMAP/"
#ml = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/25K_norm/" # Normal
ml = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/1_Bin_Std_Matrices/25K_norm/" # Binned
tempLoc = "/scratch/egadegu/ABT/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/31_Average_Diff_Interactions_AB/1_Tables/bin/"
abTreatList = ["STAG1+", "STAG2+"]
abContrList = ["STAG1-", "STAG2-"]
tadTreatList = ["STAG1+AUX", "STAG2+AUX"]
tadContrList = ["STAG1-AUX", "STAG2-AUX"]
matTreatList = ["69_127plus", "3B9_5plus"]
matContrList = ["69_127-", "3B9_5-"]

# Open File
inputFileName = fl + "1_abtb.txt"
inFile = open(inputFileName, "w")

# Loop
for i in range(0,len(abTreatList)):

  # Parameters
  abTreat = abTreatList[i]
  abContr = abContrList[i]
  tadTreat = tadTreatList[i]
  tadContr = tadContrList[i]
  matTreat = matTreatList[i]
  matContr = matContrList[i]

  # Input
  resolution = "25000"
  abTreatFileNameList = ",".join([al + abTreat + "/" + e + "_2500.bed" for e in chrList])
  abControlFileNameList = ",".join([al + abContr + "/" + e + "_2500.bed" for e in chrList])
  #treatMatrixFileName = ml + matTreat + ".txt" # Normal
  #controlMatrixFileName = ml + matContr + ".txt" # Normal
  treatMatrixFileName = ml + matTreat + "_bin.txt" # Binned
  controlMatrixFileName = ml + matContr + "_bin.txt" # Binned
  treatTadFileName = tl + tadTreat + ".bam"
  controlTadFileName = tl + tadContr + ".bam"
  temporaryLocation = tempLoc + abTreat + "/"
  outputIntraFileName = ol + abTreat + "_intra.txt"
  outputInterFileName = ol + abTreat + "_inter.txt"

  # Execution
  inFile.write(" ".join([resolution, abTreatFileNameList, abControlFileNameList, treatMatrixFileName, controlMatrixFileName, treatTadFileName, controlTadFileName, temporaryLocation, outputIntraFileName, outputInterFileName])+"\n")

inFile.close()


