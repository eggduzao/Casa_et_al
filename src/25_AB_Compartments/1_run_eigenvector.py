
# Import
import os
import sys

# Condition List
counter = 1
fl = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/25_AB_Compartments/input_eig/"
ml = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_matrix_files/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/25_AB_Compartments/1_eigenvector/"
condList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]
condNameList = ["STAG2-", "STAG2+", "STAG1-", "STAG1+"]

# Condition Loop
for i in range(0,len(condList)):

  cond = condList[i]
  condName = condNameList[i]

  # Chromosome List
  chromList = [str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Resolution List
    resList = ["25000", "50000", "100000"]

    # Resolution Loop
    for res in resList:

      # Name
      outLoc = ol+condName+"/"
      os.system("mkdir -p "+outLoc)
      name = "chr"+chrom+"_"+str(int(res)/1000)

      # Input
      normType = "KR"
      inputFileName = ml+cond+"_30.hic"
      chromN = chrom
      frag = "BP"
      resolution = res
      outputFileName = outLoc+name+".txt"

      # Write
      inFile = open(fl+str(counter)+".txt", "w")
      inFile.write("\n".join([normType, inputFileName, chromN, frag, resolution, outputFileName]))
      inFile.close()
      counter += 1


