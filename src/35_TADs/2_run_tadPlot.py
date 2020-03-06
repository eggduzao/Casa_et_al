
# Import
import os
import sys
from glob import glob

# Experiment List
jump = 10000000
step = 5000000
nameRes = 25000
res = 25000
ml = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/20_sparseToMatrix/10_5/25K_norm/"
tl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/1_TADs_GMAP/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/2_TADs_GMAP_plots/"
expList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]
nameList = ["STAG2", "STAG2", "STAG1", "STAG1"]
condList = ["-AUX", "+AUX", "-AUX", "+AUX"]

# Experiment Loop
for i in range(0,len(expList)):

  exp = expList[i]
  namey = nameList[i]
  cond = condList[i]

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

      # Input List
      inputTadList = glob(tl+exp+"/"+chrom+"/"+"*_htad.txt")

      # Input Loop
      for inputTad in inputTadList:

        # Creating output
        outLoc = ol+exp+"/"+chrom+"/"
        os.system("mkdir -p "+outLoc)

        # Name
        matName = str(s/nameRes)+"_"+str(min(s+jump,chromSizesDict[chrom])/nameRes)
        name = matName+"_"+"_".join(inputTad.split("/")[-1].split("_")[:-1])

        # Input
        start = str(s)
        end = str(min(s+jump,chromSizesDict[chrom]))
        resolution = str(res)
        namey = namey
        cond = cond
        matrixFileName = ml+exp+"/"+chrom+"/"+matName+".mat"
        tadFileName = inputTad
        newTadFileName = inputTad+"TEMP.txt"
        outputFileName = outLoc+"/"+name+".pdf"

        # Converting TADs
        print exp, chrom, name
        command = "python 2_tadPlot.py "+" ".join([start, end, resolution, namey, cond, matrixFileName, tadFileName, newTadFileName, outputFileName])
        os.system(command)


