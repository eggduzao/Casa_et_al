
# Import
import os
import sys

# Experiment List

#jump = 20000000
#step = 10000000
#nameRes = 1000000
#ml = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/20_sparseToMatrix/20_10/10K_norm/"

jump = 10000000
step = 5000000
nameRes = 1000000
ml = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/20_sparseToMatrix/10_5/10K_norm/"

#jump = 6000000
#step = 3000000
#nameRes = 1000000
#ml = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/20_sparseToMatrix/6_3/10K_norm/"

#jump = 5000000
#step = 2500000
#nameRes = 10000
#ml = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/20_sparseToMatrix/5_2.5/10K_norm/"

#jump = 5000000
#step = 5000000
#nameRes = 10000
#ml = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/20_sparseToMatrix/5_0/"

res = 10000
tl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/18_TADs/1_TADs_GMAP/10K_norm/"
dtl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/21_differential_TAD/10K_norm/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/19_TAD_plotter/10K_norm/"
expPairList = [["3B9_5-", "3B9_5plus"], ["69_127-", "69_127plus"]]
nameList = [["STAG2", "STAG2"], ["STAG1", "STAG1"]]
condList = [["-AUX", "+AUX"], ["-AUX", "+AUX"]]
expPairList = [["3B9_5-", "3B9_5plus"]]
nameList = [["STAG2", "STAG2"]]
condList = [["-AUX", "+AUX"]]

# Experiment Loop
for i in range(0,len(expPairList)):

  expPair = expPairList[i]
  namePair = nameList[i]
  condPair = condList[i]
  exp = expPair[0][:-1]

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]
  chrList = ["chr1"]

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
    startList = [40000000]

    # Start Loop
    for s in startList:

      if(s >= chromSizesDict[chrom]): break

      # Creating output
      outLoc = ol+exp+"/"+chrom+"/"
      os.system("mkdir -p "+outLoc)

      # Name
      name = "_".join([str(s/nameRes),str(min((s+jump)/nameRes,chromSizesDict[chrom]/nameRes))])

      """
      # Input
      start = str(s)
      end = str(min(s+jump,chromSizesDict[chrom]))
      resolution = str(res)
      matrixFileName = ml+expPair[0]+"/"+chrom+"/"+name+".mat"
      tadFileName = tl+expPair[0]+"/"+chrom+"_htad.txt"
      newTadFileName = tl+expPair[0]+"/"+chrom+"_htad_TEMP1.txt"
      outputFileName = ol+exp+"/"+chrom+"/"+name+".pdf"

      # Converting TADs
      command = "python 1_format_gmap.py "+" ".join([start, end, tadFileName, newTadFileName])
      os.system(command)

      # Printing matrix
      command = "R CMD BATCH '--args '"+resolution+"' '"+start+"' '"+matrixFileName+"' '"+newTadFileName+"' '"+outputFileName+" 1_tadPlotter.R 1_tadPlotter.Rout"
      os.system(command)

      # Removing TEMP files
      command = "rm "+" ".join([newTadFileName])
      os.system(command)
      """
   
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
      diffTadFileName1 = dtl+exp+"/"+chrom+"_htad_1.txt"
      diffTadFileName2 = dtl+exp+"/"+chrom+"_htad_1.txt"
      newTadFileName1 = tl+expPair[0]+"/"+chrom+"_htad_TEMP1.txt"
      newTadFileName2 = tl+expPair[1]+"/"+chrom+"_htad_TEMP2.txt"
      newDiffTadFileName1 = dtl+exp+"/"+chrom+"_dhtad_TEMP1.txt"
      newDiffTadFileName2 = dtl+exp+"/"+chrom+"_dhtad_TEMP2.txt"
      outputFileName = ol+exp+"/"+chrom+"/"+name

      # Converting TADs
      command = "python 1_format_gmap.py "+" ".join([start, end, tadFileName1, newTadFileName1])
      os.system(command)
      command = "python 1_format_gmap.py "+" ".join([start, end, tadFileName2, newTadFileName2])
      os.system(command)
      command = "python 1_format_gmap.py "+" ".join([start, end, diffTadFileName1, newDiffTadFileName1])
      os.system(command)
      command = "python 1_format_gmap.py "+" ".join([start, end, diffTadFileName2, newDiffTadFileName2])
      os.system(command)

      # Printing matrix
      command = "R CMD BATCH '--args '"+resolution+"' '"+start+"' '"+name1+"' '"+cond1+"' '"+name2+"' '"+cond2+"' '"+matrixFileName1+"' '"+matrixFileName2+"' '"+newTadFileName1+"' '"+newTadFileName2+"' '"+newDiffTadFileName1+"' '"+newDiffTadFileName2+"' '"+outputFileName+" 2_tadPlotter2.R 2_tadPlotter2.Rout"
      os.system(command)

      # Merging TADs
      command = "montage "+outputFileName+"_1.pdf "+outputFileName+"_2.pdf -tile 2x1 -geometry +0+0 "+outputFileName+".pdf"
      os.system(command)

      # Removing TEMP files
      command = "rm "+" ".join([newTadFileName1, newTadFileName2, newDiffTadFileName1, newDiffTadFileName2, outputFileName+"_1.pdf", outputFileName+"_2.pdf"])
      os.system(command)
    

