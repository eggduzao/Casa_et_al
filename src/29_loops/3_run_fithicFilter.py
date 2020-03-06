
# Import
import os
import sys

# Input List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/29_loops/2_contact_files_raw/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/29_loops/3_contact_files_filtered/"
inList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

# Input Loop
for inName in inList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # qValue List (Fetching the TOP 1 percentile)
    if(inName == "3B9_5-"): qList = ["8.385807e-28"]
    elif(inName == "3B9_5plus"): qList = ["2.0062112e-23"]
    elif(inName == "69_127-"): qList = ["1.129803e-16"] 
    elif(inName == "69_127plus"): qList = ["2.333198e-27"] 

    # qValue Loop
    for qValue in qList:

      # Parameter
      outLoc = ol + inName + "/"
      os.system("mkdir -p "+outLoc)

      # Input
      qValueThreshold = qValue
      inputContactFileName = il + inName + "/" + chrom + ".spline_pass1.res25000.significances.txt.gz"
      tempLoc = outLoc + "TEMP/"
      outputContactFileName = outLoc + chrom + ".txt"

      # Execution
      print inName, chrom
      command = "python 3_fithicFilter.py "+" ".join([qValueThreshold, inputContactFileName, tempLoc, outputContactFileName])
      os.system(command)


