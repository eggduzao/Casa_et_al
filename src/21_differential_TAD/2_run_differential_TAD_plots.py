
# Import
import os
import sys

# Experiment List
jump = 10000000
step = 5000000
nameRes = 1000000
res = 10000
ml = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/20_sparseToMatrix/10_5/10K_norm/"
tl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/18_TADs/1_TADs_GMAP/10K_norm/"
dtl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/21_differential_TAD/1_differential_TAD/10K_norm/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/21_differential_TAD/2_differential_TAD_plots/10K_norm/"
expPairList = [["3B9_5-", "3B9_5plus"], ["69_127-", "69_127plus"]]
nameList = [["STAG2", "STAG2"], ["STAG1", "STAG1"]]
condList = [["-AUX", "+AUX"], ["-AUX", "+AUX"]]

# Experiment Loop
for i in range(0,len(expPairList)):

  expPair = expPairList[i]
  namePair = nameList[i]
  condPair = condList[i]
  exp = expPair[0][:-1]

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
      name1 = namePair[0]
      cond1 = condPair[0]
      name2 = namePair[1]
      cond2 = condPair[1]
      matrixFileName1 = ml+expPair[0]+"/"+chrom+"/"+name+".mat"
      matrixFileName2 = ml+expPair[1]+"/"+chrom+"/"+name+".mat"
      tadFileName1 = tl+expPair[0]+"/"+chrom+"_htad.txt"
      tadFileName2 = tl+expPair[1]+"/"+chrom+"_htad.txt"
      diffTadFileName1 = dtl+exp+"/"+chrom+"_htad_2.txt"
      diffTadFileName2 = dtl+exp+"/"+chrom+"_htad_1.txt"
      newTadFileName1 = tl+expPair[0]+"/"+chrom+"_htad_TEMP1.txt"
      newTadFileName2 = tl+expPair[1]+"/"+chrom+"_htad_TEMP2.txt"
      newDiffTadFileName1 = dtl+exp+"/"+chrom+"_dhtad_TEMP1.txt"
      newDiffTadFileName2 = dtl+exp+"/"+chrom+"_dhtad_TEMP2.txt"
      outputFileName = ol+exp+"/"+chrom+"/"+name

      # Converting TADs
      command = "python 2_differential_TAD_plots.py "+" ".join([start, end, resolution, name1, cond1, name2, cond2, matrixFileName1, matrixFileName2, tadFileName1, tadFileName2, diffTadFileName1, diffTadFileName2, newTadFileName1, newTadFileName2, newDiffTadFileName1, newDiffTadFileName2, outputFileName])
      os.system(command)


