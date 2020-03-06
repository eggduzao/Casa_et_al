
# Import
import os
import sys

# Hic File List
chromSizesFile = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/19_Process_Losada_Hic/2_Merged_Hic/"
tl = "/scratch/egadegu/HTS/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/1_Sparse_Matrices/"
hicFileNameList = ["LCONT", "SISA1", "SISA2"]

# Open file
inputFileName = fl + "1_hts.txt"
inFile = open(inputFileName, "w")

# Hic Loop
for hicName in hicFileNameList:
  
  # Resolution List
  resList = ["25000", "50000", "100000", "250000"]

  # Resolution Loop
  for res in resList:

    # Norm List
    normList = ["none", "KR"]

    # Norm Loop
    for norm in normList:

      # Parameters
      if(res == "25000"): rName = "25K"
      elif(res == "50000"): rName = "50K"
      elif(res == "100000"): rName = "100K"
      elif(res == "250000"): rName = "250K"
      if(norm == "none"): nName = "none"
      elif(norm == "KR"): nName = "norm"
      outName = rName+"_"+nName

      # Input
      juicerCommand = "juicertools"
      kindOfMatrix = "observed"
      kindOfNormalization = norm
      unitOfResolution = "BP"
      resolution = res
      chromSizesFileName = chromSizesFile
      inputHicFileName = il + hicName + "/inter_30.hic"
      tempLocation = tl + hicName + "_" + outName + "/"
      outputMatrixFileName = ol + outName + "/" + hicName + ".txt"
      outputScoreFileName = ol + outName + "/" + hicName + "_percentiles.txt"

      # Execution
      inFile.write(" ".join([juicerCommand, kindOfMatrix, kindOfNormalization, unitOfResolution, resolution, chromSizesFileName, inputHicFileName, tempLocation, outputMatrixFileName, outputScoreFileName])+"\n")

inFile.close()


