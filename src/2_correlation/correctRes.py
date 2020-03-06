
# Import
import os
import sys
from glob import glob

# Execution
for inFileName in glob("/home/egg/Projects/hic_corr/result_all/*.txt"):
  inFileName2 = inFileName+"TEMP"
  inFile = open(inFileName,"r")
  outFile = open(inFileName2,"w")
  for line in inFile:
    ll = line.strip().split()
    if(ll[1] == "NA" or ll[1] == "NaN" or "chrX" in ll[0]): continue
    outFile.write(line)
  inFile.close()
  outFile.close()
  os.system("mv "+inFileName2+" "+inFileName)


