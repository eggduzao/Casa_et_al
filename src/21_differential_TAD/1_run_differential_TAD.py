
# Import
import os
import sys

# Experiment List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/18_TADs/1_TADs_GMAP/10K_norm/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/21_differential_TAD/1_differential_TAD/10K_norm/"
expPairList = [["3B9_5-", "3B9_5plus"], ["69_127-", "69_127plus"]]

# Experiment Loop
for expPair in expPairList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chrList:

    name = expPair[0][:-1]

    # Input
    resolution = "10000"
    maxDiffBin = "10" # This always means (x*2)+1
    tadFileName1 = il+expPair[0]+"/"+chrom+"_htad.txt"
    tadFileName2 = il+expPair[1]+"/"+chrom+"_htad.txt"
    tempLocation = ol+"TEMP/"
    outputFileName1 = ol+name+"/"+chrom+"_htad_1.txt"
    outputFileName2 = ol+name+"/"+chrom+"_htad_2.txt"

    # Converting TAD
    command = "python 1_differential_TAD.py "+" ".join([resolution, maxDiffBin, tadFileName1, tadFileName2, tempLocation, outputFileName1, outputFileName2])
    os.system(command)


