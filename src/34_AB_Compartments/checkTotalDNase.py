
# Import
import os
import sys

il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/25_AB_Compartments/2_compartments/"

# Resolution List
resList = ["50000", "100000"]

# Resolution Loop
for res in resList:

  totalDNase = 0

  # Cond List
  condList = ["STAG1-", "STAG1+", "STAG2-", "STAG2+"]

  # Condition Loop
  for cond in condList:

    # Chromosome List
    chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

    # Chromosome Loop
    for chrom in chromList:

      # Name
      inLoc = il+cond+"/"
      name = chrom+"_"+str(int(res)/1000)

      # Input
      outCountFileName = inLoc+name+".txt"
      outCountFile = open(outCountFileName, "rU")
      outCountFile.readline()
      ll = outCountFile.readline().strip().split("\t")
      totalDNase += (max(float(ll[0]),float(ll[1])) - min(float(ll[0]),float(ll[1])))
      outCountFile.close()     

  print res, totalDNase


