
# Import
import os
import sys

# Peak List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/expression/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/2_active_inactive_genes/"
peakList = ["STAG1", "STAG2"]

# Peak Loop
for peakName in peakList:

  # Input
  expThresh = "20"
  expressionFileName = il + peakName + "_minus.txt"
  outputActiveFileName = ol + peakName + "_active.txt"
  outputInactiveFileName = ol + peakName + "_inactive.txt"

  # Execution
  command = "python 2_listOfActiveInactiveGenes.py "+" ".join([expThresh, expressionFileName, outputActiveFileName, outputInactiveFileName])
  os.system(command)


