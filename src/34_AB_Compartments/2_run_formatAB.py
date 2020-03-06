
# Import
import os
import sys

# Condition List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/3_AB_Compartments/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/3_AB_Compartments/1_eigenvector/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/3_AB_Compartments/2_compartments/"
condList = ["STAG1-", "STAG1+", "STAG2-", "STAG2+"]

# Write
inFileName = fl + "2_fab.txt"
inFile = open(inFileName, "w")

# Condition Loop
for cond in condList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Resolution List
    #resList = ["50000", "100000"]
    resList = ["1000000", "2500000"]

    # Resolution Loop
    for res in resList:

      # Name
      inLoc = il+cond+"/"
      outLoc = ol+cond+"/"
      os.system("mkdir -p "+outLoc)
      name = chrom+"_"+str(int(res)/1000)

      # Input
      chromosome = chrom
      resolution = res
      chromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
      enhChromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.enhanced"
      eigenFileName = inLoc+name+".txt"
      dnaseFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Data/dnase/bam/HCT116_DNase-seq_UW.bam"
      tempLocation = "/scratch/egadegu/AB/" + name + "/"
      outBedFileName = outLoc+name+".bed"
      outBwFileName = outLoc+name+".bw"
      outCountFileName = outLoc+name+".txt"

      # Write
      inFile.write(" ".join([chromosome, resolution, chromSizesFileName, enhChromSizesFileName, eigenFileName, dnaseFileName, tempLocation, outBedFileName, outBwFileName, outCountFileName])+"\n")
      
inFile.close()


