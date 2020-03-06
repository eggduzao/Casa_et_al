
# Import
import os
import sys

# Experiment List
jump = 10000000
step = 5000000
nameRes = 25000
res = 25000
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_matrix_files/25K_norm/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/20_sparseToMatrix/10_5/25K_norm/"
expList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

# Experiment Loop
for exp in expList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Chromosome Dictionary
    chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
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
      outLoc = ol+exp+"/"+chrom+"/"
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
      command = "python 1_sparseToMatrix.py "+" ".join([start, end, resolution, sparseFileName, outputFileName])
      os.system(command)
    

