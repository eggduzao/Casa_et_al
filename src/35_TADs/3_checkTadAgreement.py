
# Import
import os
import sys

# Input
inputTableFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/0_total_TAD_counts/number_tads.txt"

# Execution
inputTableFile = open(inputTableFileName, "rU")
inputTableFile.readline()
for line in inputTableFile:
  ll = line.strip().split("\t")
  myList = [int(e) for e in ll[1:]]
  delta = 0
  for i in range(0, len(myList)-1):
    for j in range(i+1, len(myList)):
      delta += abs(myList[i] - myList[j])
  print ll[0], delta
inputTableFile.close()


