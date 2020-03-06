
# Import
import os
import sys

# Input List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/5_Loops/2_contact_files_raw/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/5_Loops/3_contact_files_filtered/"
tl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/5_Loops/4_percentiles/"
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
      outputFolder = ol + "_".join([inName, newLb, newUb]) + "/"
      os.system("mkdir -p "+outputFolder)
      thresholdFileName = tl + "_".join([inName, newLb, newUb]) + ".txt"
      print inName, newLb, newUb

      # Fetching p-value
      pValue = 0.01
      thresholdFile = open(thresholdFileName, "rU")
      for line in thresholdFile:
        ll = line.strip().split("\t")
        if(ll[0] == "1"):
          pValue = ll[1]
      thresholdFile.close()

      # Chromosome List
      chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

      # Chromosome Loop
      for chrom in chromList:

        # Input
        pValueThreshold = pValue
        percentileThreshold = "80"
        inputContactFileName = il + inputFolder + "/" + chrom + ".txt.gz"
        tempLoc = outputFolder + "TEMP/"
        outputContactFileName = outputFolder + chrom + ".txt"

        # Execution
        command = "python 3_fithicFilter.py "+" ".join([pValueThreshold, percentileThreshold, inputContactFileName, tempLoc, outputContactFileName])
        os.system(command)

      # Contatenating chromosomes
      tempCat = "./tempCat.bed"
      command = "cat " + outputFolder + "chr*.txt > " + tempCat
      os.system(command)

      # Sorted loops
      finalFileName = ol + "_".join([inName, newLb, newUb]) + ".txt"
      command = "sort -k1,1 -k2,2n " + tempCat + " > " + finalFileName
      os.system(command)

      # Removing temporary file
      command = "rm "+tempCat
      os.system(command)


