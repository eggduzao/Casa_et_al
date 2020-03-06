
# Import
import os
import sys

# Condition List
counter = 1
fl = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/25_AB_Compartments/input_ab/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/25_AB_Compartments/1_eigenvector/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/25_AB_Compartments/2_compartments/"
condList = ["STAG1-", "STAG1+", "STAG2-", "STAG2+"]

# Condition Loop
for cond in condList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Resolution List
    resList = ["50000", "100000"]

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
      chromSizesFileName = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
      enhChromSizesFileName = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.enhanced"
      eigenFileName = inLoc+name+".txt"
      dnaseFileName = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/dnase/bam/HCT116_DNase-seq_UW.bam"
      tempLocation = "/scratch/eduardo/AB"+str(counter)+"/"
      outBedFileName = outLoc+name+".bed"
      outBwFileName = outLoc+name+".bw"
      outCountFileName = outLoc+name+".txt"

      # Write
      inFile = open(fl+str(counter)+".txt", "w")
      inFile.write("\n".join([chromosome, resolution, chromSizesFileName, enhChromSizesFileName, eigenFileName, dnaseFileName, tempLocation, outBedFileName, outBwFileName, outCountFileName]))
      inFile.close()
      counter += 1


