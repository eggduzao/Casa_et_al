
# Import
import os
import sys

# Experiment List
prefix = "T_3_80_25_100_5_10_0.95_0.5"
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/1_TADs_GMAP/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/26_differential_TAD_25K/1_differential_TAD/"
expPairList = [["3B9_5-", "3B9_5plus"], ["69_127-", "69_127plus"], ["69_127-", "3B9_5-"]]

# Experiment Loop
for expPair in expPairList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chrList:

    if(expPair[1] != "3B9_5-"): name = expPair[0][:-1]
    else: name = "WT"

    # Input
    resolution = "25000"
    maxDiffBin = "4" # This always means (x*2)+1
    tadFileName1 = il+expPair[0]+"/"+chrom+"/"+prefix+"_tad.txt"
    tadFileName2 = il+expPair[1]+"/"+chrom+"/"+prefix+"_tad.txt"
    tempLocation = ol+"TEMP/"
    outputFileName1 = ol+name+"/"+chrom+"_tad_1.txt"
    outputFileName2 = ol+name+"/"+chrom+"_tad_2.txt"

    # Converting TAD
    command = "python 1_differential_TAD.py "+" ".join([resolution, maxDiffBin, tadFileName1, tadFileName2, tempLocation, outputFileName1, outputFileName2])
    os.system(command)


