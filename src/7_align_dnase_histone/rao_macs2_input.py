# Import
import os
import sys

# Bam list
counter = 1
fl = "./input_macs_rao/"
contrFileName = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/rao/info/experiments_and_controls.txt"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/rao/bam/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/rao/macs/"
treatList = ["CHIP001_RAD21_untreated_SE", "CHIP002_RAD21_treated_SE", "CHIP003_SMC1_untreated_PE", "CHIP004_SMC1_treated_PE", "CHIP005_CTCF_untreated_PE", "CHIP006_CTCF_treated_PE", "CHIP007_CTCF_untreated_SE", "CHIP008_CTCF_treated_SE", "CHIP009_H3K27Ac_untreated_SE", "CHIP010_H3K27Ac_treated_SE", "CHIP011_H3K4me3_untreated_SE", "CHIP012_H3K4me3_treated_SE", "CHIP013_H3K4me1_untreated_SE", "CHIP014_H3K4me1_treated_SE", "CHIP015_H3K36me3_untreated_SE", "CHIP016_H3K36me3_treated_SE", "CHIP017_H3K27me3_untreated_SE", "CHIP018_H3K27me3_treated_SE", "CHIP019_H3K9me3_untreated_SE", "CHIP020_H3K9me3_treated_SE", "CHIP021_H4K16Ac_untreated_SE", "CHIP022_H4K16Ac_treated_SE", "CHIP023_H3K79me2_untreated_SE", "CHIP024_H3K79me2_treated_SE", "CHIP025_H4K20me3_untreated_SE", "CHIP026_H4K20me3_treated_SE", "CHIP027_H2AZ_untreated_SE", "CHIP028_H2AZ_treated_SE", "CHIP029_NIPBL_untreated_SE", "CHIP030_NIPBL_treated_SE"]

# Control dictionary
contrDict = dict()
controlFile = open(contrFileName, "rU")
for line in controlFile:
  ll = line.strip().split("\t")
  contrDict[ll[0]] = ll[1]
controlFile.close()

# Fastq Loop
for treat in treatList:

  # Auxiliary
  if("_PE" in treat): dtype = "PE"
  else: dtype = "ChIP-seq"

  # Fetching control
  contr = contrDict[treat.split("_")[0]]

  # Parameters
  dataType = dtype
  treatmentFileName = il+treat+".bam"
  controlFileName = il+contr+".bam"
  tempLocation = "/scratch/eduardo/stag_macs_rao/" + treat + "/"
  outputLocation = ol

  # Creating files
  inFileName = fl+str(counter)+".txt"
  inFile = open(inFileName,"w")
  inFile.write("\n".join([dataType, treatmentFileName, controlFileName, tempLocation, outputLocation]))
  inFile.close()
  counter += 1


