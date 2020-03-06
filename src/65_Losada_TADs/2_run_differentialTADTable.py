
# Import
import os
import sys
from glob import glob

# Experiment List
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/36_Losada_TADs/1_TADs_GMAP/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/36_Losada_TADs/2_Diff_TAD_Table/"
expPairList = [["SISA1", "LCONT"], ["SISA2", "LCONT"], ["SISA1", "SISA2"]]

# Experiment Loop
for expPair in expPairList:

  # Threshold List
  threshList = ["0", "25000", "50000", "75000", "100000", "125000"]
  threshLabel = ["0_bins", "1_bins", "2_bins", "3_bins", "4_bins", "5_bins"]

  # Thresold Loop
  for i in range(0,len(threshList)):

    # Chromosome List
    chrList = ["chr"+str(e) for e in range(1,23)+["X"]]
    name = expPair[0] + "_" + expPair[1]

    # Parameters
    threshNb = threshList[i]
    threshName = threshLabel[i]
    outLoc = ol + threshName + "/"
    os.system("mkdir -p " + outLoc)

    # Chromosome Loop
    for chrom in chrList:

      # Input
      maxDiffInBp = threshNb
      tadFileName1 = il + expPair[0] + "_" + chrom + "_tad.txt"
      tadFileName2 = il + expPair[1] + "_" + chrom + "_tad.txt"
      tempLocation = outLoc + chrom + "/"
      outputFileName = outLoc + chrom + ".txt"

      # Execution
      command = "python /usr/users/egadegu/Projects/Wendt_Stag/Code/36_Losada_TADs/2_differentialTADTable.py "+" ".join([maxDiffInBp, tadFileName1, tadFileName2, tempLocation, outputFileName])
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
    command = "rm -rf "+catFileName
    os.system(command)


