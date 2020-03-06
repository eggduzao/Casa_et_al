
# Import
import os
import sys
from glob import glob

# Experiment List
prefix = "T_3_80_25_100_5_10_0.95_0.5"
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/1_TADs_GMAP/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/27_differential_TAD_table/"
expPairList = [["3B9_5-", "3B9_5plus"], ["69_127-", "69_127plus"], ["69_127-", "3B9_5-"]]

# Experiment Loop
for expPair in expPairList:

  # Threshold List
  threshList = ["0", "25000", "50000", "75000", "100000", "125000"]
  threshLabel = ["0_bins", "1_bins", "2_bins", "3_bins", "4_bins", "5_bins"]

  # Thresold Loop
  for i in range(0,len(threshList)):

    # Chromosome List
    chrList = ["chr"+str(e) for e in range(1,23)+["X"]]
    if(expPair[1] != "3B9_5-"): name = expPair[0][:-1]
    else: name = "WT"

    # Parameters
    threshNb = threshList[i]
    threshName = threshLabel[i]
    outLoc = ol + threshName + "/" + name + "/"

    # Chromosome Loop
    for chrom in chrList:

      # Input
      maxDiffInBp = threshNb
      tadFileName1 = il + expPair[0] + "/" + chrom + "/" + prefix + "_tad.txt"
      tadFileName2 = il + expPair[1] + "/" + chrom + "/" + prefix + "_tad.txt"
      tempLocation = outLoc + chrom + "/"
      outputFileName = outLoc + chrom + ".txt"

      # Execution
      command = "python 1_differential_TAD.py "+" ".join([maxDiffInBp, tadFileName1, tadFileName2, tempLocation, outputFileName])
      os.system(command)

    # Cat
    catFileName = ol + threshName + "/" + name + "_cat.txt"
    command = "cat "+" ".join(glob(outLoc + "chr*.txt"))+" > "+ catFileName
    os.system(command)
 
    # Sort
    sortFileName = ol + threshName + "/" + name + ".txt"
    command = "sort -k1,1 -k2,2n "+catFileName+" > "+sortFileName
    os.system(command)

    # Insert header
    header = "\t".join(["CHROMOSOME", "TAD1_START", "TAD1_END", "TAD2_START", "TAD2_END", "STATE_NAME", "STATE_NUMBER"])+"\n"
    command = "sed -i '1 i\\"+header+"' "+sortFileName
    os.system(command)

    # Remove files
    command = "rm -rf "+outLoc+" "+catFileName
    os.system(command)


