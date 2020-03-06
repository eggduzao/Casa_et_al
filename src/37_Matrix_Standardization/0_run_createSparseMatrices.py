
## Import
import os
import sys

# Hic file list
il = "/home/egg/Projects/Wendt_Stag/Data/stag_matrix_files/"
ol = "/home/egg/Projects/Wendt_Stag/Data/stag_matrix_files/"
hicNameList = ["3B9_5-_30", "3B9_5plus_30", "69_127-_30", "69_127plus_30"]

# Hic Loop
for hicName in hicNameList:
  
  # Resolution List
  resList = ["100000", "250000", "500000", "1000000"]
  resList = ["250000", "500000", "1000000"]

  # Resolution Loop
  for res in resList:

    # Norm List
    normList = ["none", "KR"]

    for norm in normList:

      # Parameters
      if(res == "10000"): rName = "10K"
      elif(res == "25000"): rName = "25K"
      elif(res == "50000"): rName = "50K"
      elif(res == "100000"): rName = "100K"
      elif(res == "250000"): rName = "250K"
      elif(res == "500000"): rName = "500K"
      elif(res == "1000000"): rName = "1000K"
      if(norm == "none"): nName = "none"
      elif(norm == "KR"): nName = "norm"
      outName = rName+"_"+nName
      hName = "_".join(hicName.split("_")[:-1])
      outLoc = ol + outName + "/"
      os.system("mkdir -p "+outLoc)

      # Input
      juicerCommand = "juicertools"
      kindOfMatrix = "observed"
      kindOfNormalization = norm
      unitOfResolution = "BP"
      resolution = res
      chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19"
      inputHicFileName = il + hicName + ".hic"
      tempLocation = "./TEMP/"
      outputMatrixFileName = outLoc + hName + ".txt"
      outputScoreFileName = outLoc + hName + "_percentiles.txt"

      # Execution
      command = "python 0_createSparseMatrices.py "+" ".join([juicerCommand, kindOfMatrix, kindOfNormalization, unitOfResolution, resolution, chromSizesFileName, inputHicFileName, tempLocation, outputMatrixFileName, outputScoreFileName])
      os.system(command)


