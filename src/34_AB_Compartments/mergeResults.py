
# Import
import os
import sys

# Input
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/25_AB_Compartments/2_compartments/"
enhChromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.enhanced"

# Condition List
condList = ["STAG1-", "STAG1+", "STAG2-", "STAG2+"]

# Condition Loop
for cond in condList:

  # Output files
  outputBwFileName = il+cond+".bw"
  tempWigFileName = il+cond+".wig"
  tempWigFile = open(tempWigFileName, "w")
  outputBedFileName = il+cond+".bed"
  outputBedFile = open(outputBedFileName, "w")

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # Resolution List
    resList = ["50000"]

    # Resolution Loop
    for res in resList:

      # Name
      inLoc = il+cond+"/"
      name = chrom+"_"+str(int(res)/1000)

      # Input
      inBedFileName = inLoc+name+".bed"
      inBwFileName = inLoc+name+".bw"

      # Fetching bed
      inBedFile = open(inBedFileName, "rU")
      for line in inBedFile: outputBedFile.write(line)
      inBedFile.close()

      # Fetching wig
      zomgWigFileName = il+cond+"_ZOMG.wig"
      command = "bigWigToWig "+" ".join([inBwFileName, zomgWigFileName])
      os.system(command)
      zomgWigFile = open(zomgWigFileName, "rU")
      for line in zomgWigFile: tempWigFile.write(line)
      zomgWigFile.close()
      os.system("rm "+zomgWigFileName)

  # Close
  tempWigFile.close()
  outputBedFile.close()
  command = "wigToBigWig "+" ".join([tempWigFileName, enhChromSizesFileName, outputBwFileName])
  os.system(command)
  os.system("rm "+tempWigFileName)


