
# Import
import os
import sys

# Condition List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/input/"
ml1 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/rao_new_untreated_hic/"
ml2 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/41_FigS9A_AB_RaoNeg/1_eigenvector/"
#condList = ["Untreated", "Untreated_Synchronized", "3B9_5-", "69_127-"]
#condNameList = ["Untreated", "Untreated_Synchronized", "STAG2-", "STAG1-"]

condList = ["3B9_5plus", "69_127plus"]
condNameList = ["STAG2+", "STAG1+"]

# Write
inFileName = fl + "1_eig.txt"
inFile = open(inFileName, "w")

# Condition Loop
for i in range(0,len(condList)):

  cond = condList[i]
  condName = condNameList[i]

  # Chromosome List
  chromList = [str(e) for e in list(range(1,23))+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Resolution List
    resList = ["100000"]

    # Resolution Loop
    for res in resList:

      # Name
      outLoc = ol + condName + "/"
      os.system("mkdir -p " + outLoc)
      name = "chr" + chrom

      # Parameters
      if("treated" in cond):
        ml = ml1
        suff = "."
      else:
        ml = ml2
        suff = "_30."

      # Input
      normType = "KR"
      inputFileName = ml + cond + suff + "hic"
      chromN = chrom
      frag = "BP"
      resolution = res
      outputFileName = outLoc + name + ".txt"

      # Write
      inFile.write(" ".join([normType, inputFileName, chromN, frag, resolution, outputFileName])+"\n")
      
inFile.close()


