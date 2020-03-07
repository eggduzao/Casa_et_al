
# Import
import os
import sys

# Condition List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/41_FigS9A_AB_RaoNeg/2_compartments/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/41_FigS9A_AB_RaoNeg/3_AB_Plot/"
bloomFolderList = ["STAG1+", "STAG2+", "STAG1+", "STAG2+", "STAG1-", "STAG2-", "STAG1-", "STAG2-"]
contrFolderList = ["Untreated", "Untreated", "Untreated_Synchronized", "Untreated_Synchronized", "Untreated", "Untreated", "Untreated_Synchronized", "Untreated_Synchronized"]
outNameList = ["STAG1+_v_Untreated", "STAG2+_v_Untreated", "STAG1+_v_Synchronized", "STAG2+_v_Synchronized", "STAG1-_v_Untreated", "STAG2-_v_Untreated", "STAG1-_v_Synchronized", "STAG2-_v_Synchronized"]
xLabelList = ["STAG1+", "STAG2+", "STAG1+", "STAG2+", "STAG1-", "STAG2-", "STAG1-", "STAG2-"]
yLabelList = ["Untreated", "Untreated", "Synchronized", "Synchronized", "Untreated", "Untreated", "Synchronized", "Synchronized"]

# Write
inFileName1 = fl + "31_cpt.txt"
inFile1 = open(inFileName1, "w")
inFileName2 = fl + "32_cpt.txt"
inFile2 = open(inFileName2, "w")

# Condition Loop
for i in range(0,len(bloomFolderList)):

  # Input File List
  bloom_folder = bloomFolderList[i]
  contr_folder = contrFolderList[i]
  out_name = outNameList[i]

  # Parameters
  chrList = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
             "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
  command = "cat " + " ".join([il + bloom_folder + "/" + e + ".bed" for e in chrList]) + " > " + il + bloom_folder + "/" + out_name + ".bed"
  os.system(command)
  command = "cat " + " ".join([il + contr_folder + "/" + e + ".bed" for e in chrList]) + " > " + il + contr_folder + "/" + out_name + ".bed"
  os.system(command)

  # Input 1
  abTreatFileName = il + bloom_folder + "/" + out_name + ".bed"
  abContrFileName = il + contr_folder + "/" + out_name + ".bed"
  outputFileName = il + out_name + ".bed"
  outputBedFilePrefix = il + out_name

  # Write 1
  inFile1.write(" ".join([abTreatFileName, abContrFileName, outputFileName, outputBedFilePrefix])+"\n")

  # Input 2
  xLabel = xLabelList[i]
  yLabel = yLabelList[i]
  inputFileName = il + out_name + ".bed"
  outputFileName = ol + out_name + ".pdf"

  # Write 2
  inFile2.write(" ".join([xLabel, yLabel, inputFileName, outputFileName])+"\n")

inFile1.close()
inFile2.close()


