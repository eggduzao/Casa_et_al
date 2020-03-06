
# Import
import os
import sys

# Input List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/29_loops/2_contact_files_raw/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/29_loops/4_percentiles/"
inList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

# Input Loop
for inName in inList:

  # Input
  inputLocation = il + inName + "/"
  tempLoc = ol + "TEMP/"
  outputFileName = ol + inName + ".txt"

  # Execution
  print inName
  command = "python 4_calculatingScorePercentiles.py "+" ".join([inputLocation, tempLoc, outputFileName])
  os.system(command)


