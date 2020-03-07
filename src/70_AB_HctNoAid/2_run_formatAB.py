
# Import
import os
import sys

# Condition List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/41_FigS9A_AB_RaoNeg/1_eigenvector/"
tl = "/usr/users/egadegu/scratch/FAB/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/41_FigS9A_AB_RaoNeg/2_compartments/"
#condList = ["Untreated", "Untreated_Synchronized", "STAG2-", "STAG1-"]
condList = ["STAG2+", "STAG1+"]

# Write
inFileName = fl + "2_fab.txt"
inFile = open(inFileName, "w")

# Condition Loop
for cond in condList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in list(range(1,23))+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Resolution List
    resList = ["100000"]

    # Resolution Loop
    for res in resList:

      # Name
      inLoc = il + cond + "/"
      outLoc = ol + cond + "/"
      os.system("mkdir -p "+ outLoc)
      name = chrom

      # Input
      chromosome = chrom
      resolution = res
      chromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
      enhChromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.enhanced"
      eigenFileName = inLoc + name + ".txt"
      dnaseFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Data/dnase/bam/HCT116_DNase-seq_UW.bam"
      tempLocation = tl + cond + "/" + name + "/"
      outBedFileName = outLoc + name + ".bed"
      outBwFileName = outLoc + name + ".bw"
      outCountFileName = outLoc + name + ".txt"

      # Write
      inFile.write(" ".join([chromosome, resolution, chromSizesFileName, enhChromSizesFileName, eigenFileName, dnaseFileName, tempLocation, outBedFileName, outBwFileName, outCountFileName])+"\n")
      
inFile.close()


