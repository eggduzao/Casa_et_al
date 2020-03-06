
## Import
import os
import sys

# Hic file list
hicFileNameList = ["/media/egg/sbc/AG_Papantonis/HiC-STAG1_2/merged/3B9_5-/3B9_5-_30.hic", "/media/egg/sbc/AG_Papantonis/HiC-STAG1_2/merged/3B9_5plus/3B9_5plus_30.hic", "/media/egg/sbc/AG_Papantonis/HiC-STAG1_2/merged/69_127-/69_127-_30.hic", "/media/egg/sbc/AG_Papantonis/HiC-STAG1_2/merged/69_127plus/69_127plus_30.hic"]
ol = "/media/egg/sbc/AG_Papantonis/Eduardo/Papantonis_Stag/Data/stag_matrix_files/"

# Hic Loop
for hicFileName in hicFileNameList:
  
  # Resolution List
  resList = ["10000", "25000", "50000"]

  # Resolution Loop
  for res in resList:

    # Norm List
    normList = ["none", "KR"]

    for norm in normList:

      # Parameters
      hicName = hicFileName.split("/")[-2]
      if(res == "10000"): rName = "10K"
      elif(res == "25000"): rName = "25K"
      elif(res == "50000"): rName = "50K"
      if(norm == "none"): nName = "none"
      elif(norm == "KR"): nName = "norm"
      outName = rName+"_"+nName
      outLoc = "/home/egg/Projects/Papantonis_Stag/Data/stag_matrix_files/25K_none/"

      # Input
      juicerCommand = "juicertools"
      kindOfMatrix = "observed"
      kindOfNormalization = norm
      unitOfResolution = "BP"
      resolution = res
      chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19"
      inputHicFileName = hicFileName
      tempLocation = ol+outName+"/TEMP/"
      outputMatrixFileName = ol+outName+"/"+hicName+".txt"
      outputScoreFileName = ol+outName+"/"+hicName+"_percentiles.txt"

      # Execution
      command = "python createInputMatrices.py "+" ".join([juicerCommand, kindOfMatrix, kindOfNormalization, unitOfResolution, resolution, chromSizesFileName, inputHicFileName, tempLocation, outputMatrixFileName, outputScoreFileName])
      os.system(command)


