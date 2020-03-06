
# Import
import os
import sys

# Input List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/5_Loops/2_contact_files_raw/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/5_Loops/4_percentiles/"
inList = ["STAG1_minus_auxin", "STAG1_plus_auxin", "STAG2_minus_auxin", "STAG2_plus_auxin"]

# Input Loop
for inName in inList:

  # lb List
  lbList = ["50000", "150000"]

  # lb Loop
  for lb in lbList:

    # ub List
    ubList = ["5000000", "10000000"]

    # ub Loop
    for ub in ubList:

      # Parameters
      inputFolder = "_".join([inName, lb, ub])
      newLb = lb[:-3] + "Kb"
      newUb = ub[:-6] + "Mb"
      outputName = "_".join([inName, newLb, newUb])

      # Input
      inputLocation = il + inputFolder + "/"
      tempLoc = ol + "TEMP/"
      outputFileName = ol + outputName + ".txt"

      # Execution
      print outputName
      command = "python 4_calculatingScorePercentiles.py "+" ".join([inputLocation, tempLoc, outputFileName])
      os.system(command)
      break
    break
  break

