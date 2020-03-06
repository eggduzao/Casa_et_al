
# Import
import os
import sys

# Experiment List
jump = 10000000
step = 5000000
nameRes = 25000
res = 25000
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/1_Sparse_Matrices/25K_norm/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/36_Losada_TADs/0_Sparse_To_Matrix/25K_norm/"
matrixList = ["LCONT", "SISA1", "SISA2"]

# Experiment Loop
for exp in matrixList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Chromosome Dictionary
    chromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
    chromSizesDict = dict()
    chromSizesFile = open(chromSizesFileName,"rU")
    for line in chromSizesFile:
      ll = line.strip().split("\t")
      chromSizesDict[ll[0]] = int(ll[1])
    chromSizesFile.close()
    chromList = chromSizesDict.keys()

    # Start List
    startList = range(0,chromSizesDict[chrom]+step,step)
    
    # Start Loop
    for s in startList:
 
      if(s >= chromSizesDict[chrom]): break

      # Creating output
      outLoc = ol + exp + "/" + chrom + "/"
      os.system("mkdir -p "+outLoc)

      # Name
      name = "_".join([str(s/nameRes),str(min((s+jump)/nameRes,chromSizesDict[chrom]/nameRes))])

      # Input
      start = str(s)
      end = str(min(s+jump,chromSizesDict[chrom]))
      resolution = str(res)
      sparseFileName = il+exp+"/"+chrom+".txt"
      outputFileName = ol+exp+"/"+chrom+"/"+name+".mat"

      # Execution
      print exp, chrom, name
      command = "python /usr/users/egadegu/Projects/Wendt_Stag/Code/36_Losada_TADs/0_2_sparseToMatrix.py "+" ".join([start, end, resolution, sparseFileName, outputFileName])
      os.system(command)
    

