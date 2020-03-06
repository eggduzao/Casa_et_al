
# Import
import os
import sys

# Input location
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/18_TADs/1_TADs_GMAP/10K_norm/"
outTadFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/18_TADs/1_TADs_GMAP/10K_norm/number_tads.txt"
outhTadFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/18_TADs/1_TADs_GMAP/10K_norm/number_htads.txt"
outthTadFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/18_TADs/1_TADs_GMAP/10K_norm/number_all.txt"

# Execution
nameList = ["69_127plus", "69_127-", "3B9_5plus", "3B9_5-"]
header = ["CHROM", "STAG1+", "STAG1-", "STAG2+", "STAG2-"]
ttotalCount = [0, 0, 0, 0]
htotalCount = [0, 0, 0, 0]
thtotalCount = [0, 0, 0, 0]
outTadFile = open(outTadFileName,"w")
outhTadFile = open(outhTadFileName,"w")
outthTadFile = open(outthTadFileName,"w")
outTadFile.write("\t".join(header)+"\n")
outhTadFile.write("\t".join(header)+"\n")
outthTadFile.write("\t".join(header)+"\n")
for chrom in ["chr"+str(e) for e in range(1,23)+["X"]]:
  tCount = [0, 0, 0, 0]
  hCount = [0, 0, 0, 0]
  thCount = [0, 0, 0, 0]
  for i in range(0,len(nameList)):
    inFileName = il+nameList[i]+"/"+chrom+"_htad.txt"
    inFile = open(inFileName,"rU")
    inFile.readline()
    tadCounter = 0
    htadCounter = 0
    for line in inFile:
      ll = line.strip().split("\t")
      if(ll[2] == "1"):
        tCount[i] += 1
        thCount[i] += 1
        ttotalCount[i] += 1
        thtotalCount[i] += 1
      elif(ll[2] == "2"):
        hCount[i] += 1
        thCount[i] += 1
        htotalCount[i] += 1
        thtotalCount[i] += 1
    inFile.close()
  outTadFile.write("\t".join([chrom]+[str(e) for e in tCount])+"\n")
  outhTadFile.write("\t".join([chrom]+[str(e) for e in hCount])+"\n")
  outthTadFile.write("\t".join([chrom]+[str(e) for e in thCount])+"\n")
outTadFile.write("\t".join(["TOTAL"]+[str(e) for e in ttotalCount])+"\n")
outhTadFile.write("\t".join(["TOTAL"]+[str(e) for e in htotalCount])+"\n")
outthTadFile.write("\t".join(["TOTAL"]+[str(e) for e in thtotalCount])+"\n")

# Closing files
outTadFile.close()
outhTadFile.close()
outthTadFile.close()


