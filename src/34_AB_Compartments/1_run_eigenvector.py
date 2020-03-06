
# Import
import os
import sys

# Condition List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/3_AB_Compartments/input/"
ml = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/3_AB_Compartments/1_eigenvector/"
condList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]
condNameList = ["STAG2-", "STAG2+", "STAG1-", "STAG1+"]

# Write
inFileName = fl + "1_eig.txt"
inFile = open(inFileName, "w")

# Condition Loop
for i in range(0,len(condList)):

  cond = condList[i]
  condName = condNameList[i]

  # Chromosome List
  chromList = [str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Resolution List
    #resList = ["25000", "50000", "100000"]
    resList = ["1000000", "2500000"]

    # Resolution Loop
    for res in resList:

      # Name
      outLoc = ol+condName+"/"
      os.system("mkdir -p "+outLoc)
      name = "chr"+chrom+"_"+str(int(res)/1000)

      # Input
      normType = "KR"
      inputFileName = ml + cond + "_30.hic"
      chromN = chrom
      frag = "BP"
      resolution = res
      outputFileName = outLoc + name + ".txt"

      # Write
      inFile.write(" ".join([normType, inputFileName, chromN, frag, resolution, outputFileName])+"\n")
      
inFile.close()


