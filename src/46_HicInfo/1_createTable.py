
# Import
import os
import sys
from glob import glob

# Input
il = "/home/egg/Projects/Wendt_Stag/Results/17_HicInfo/1_Input/"
tableFileName = "/home/egg/Projects/Wendt_Stag/Results/17_HicInfo/2_Table/table.tsv"
fileList = ["stag1_aux_rep1", "stag1_aux_rep2", "stag1_noaux_rep1", "stag1_noaux_rep2", "stag2_aux_rep1", "stag2_aux_rep2", "stag2_noaux_rep1", "stag2_noaux_rep2"]

# Create dictionary
tableDict = dict()
flagFirst = True
firstVec = ["Dataset"]
for inFileN in fileList:
  inFileName = il + inFileN + ".txt"
  inFile = open(inFileName, "rU")
  vec = [inFileN]
  for line in inFile:
    ll = line.strip().split(": ")
    vec.append(ll[1].replace(" ",""))
    if(flagFirst): firstVec.append("".join(["\"",ll[0],"\""]))
  inFile.close()
  tableDict[inFileN] = vec
  flagFirst = False

# Print table
tableFile = open(tableFileName, "w")
tableFile.write("\t".join(firstVec)+"\n")
for inFileN in fileList: tableFile.write("\t".join(tableDict[inFileN])+"\n")
tableFile.close()


