
# Import
import os
import sys

# Condition List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/23_AB_Compartment_Losada/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/23_AB_Compartment_Losada/1_Eigenvector/"
tl = "/scratch/egadegu/FAB/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/23_AB_Compartment_Losada/2_Compartments/"
condList = ["LCONT", "SISA1", "SISA2"]

# Opening file
inFileName = fl + "2_fab.txt"
inFile = open(inFileName,"w")

# Condition Loop
for cond in condList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Name
    inLoc = il + cond + "/"
    outLoc = ol + cond + "/"
    os.system("mkdir -p "+outLoc)

    # Input
    chromosome = chrom
    resolution = "100000"
    chromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
    enhChromSizesFileName = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.enhanced"
    eigenFileName = inLoc + chrom + ".txt"
    dnaseFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Data/losada_new/dnase/MCF10_EtOH.bam"
    tempLocation = tl + cond + "_" + chrom + "/"
    outBedFileName = outLoc + chrom + ".bed"
    outBwFileName = outLoc + chrom + ".bw"
    outCountFileName = outLoc + chrom + ".txt"

    # Write
    inFile.write(" ".join([chromosome, resolution, chromSizesFileName, enhChromSizesFileName, eigenFileName, dnaseFileName, tempLocation, outBedFileName, outBwFileName, outCountFileName])+"\n")

inFile.close()


