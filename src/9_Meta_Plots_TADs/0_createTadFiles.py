
# Import
import os
import sys

# Input
il = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/4_TADs/6_differential_TAD_tables/2_bins/"
ol = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/8_Meta_Plots/3_Tad_Files/"
diffTadTableList = ["3B9_5", "69_127", "WT"]

# Parameters
tempLoc = "./TEMP/"
command = "mkdir -p "+tempLoc
os.system(command)

# Loop
for diffName in diffTadTableList:

  # Stag name
  if(diffName == "3B9_5"): stag = "STAG2"
  elif(diffName == "69_127"): stag = "STAG1"
  else: stag = diffName

  # Creating diff tad dictionary
  diffTadDict = dict() # NUMBER_NAME -> [CHROM, TAD1_START, TAD1_END, TAD2_START, TAD2_END]
  diffTadTableFileName = il + diffName + ".txt"
  diffTadTableFile = open(diffTadTableFileName, "rU")
  diffTadTableFile.readline()
  for line in diffTadTableFile:
    ll = line.strip().split("\t")
    key = "_".join([ll[6], ll[5]])
    try: diffTadDict[key].append(ll[:5])
    except Exception: diffTadDict[key] = [ll[:5]]
  diffTadTableFile.close()

  # Writing TAD files
  for key in diffTadDict.keys():
    outputTad1FileName = ol + "_".join([key, stag, "TAD1"]) + ".txt"
    outputTad2FileName = ol + "_".join([key, stag, "TAD2"]) + ".txt"
    outputTad1File = open(outputTad1FileName, "w")
    outputTad2File = open(outputTad2FileName, "w")
    for vec in diffTadDict[key]:
      outputTad1File.write("\t".join(vec[:3])+"\n")
      outputTad2File.write("\t".join([vec[0], vec[3], vec[4]])+"\n")
    outputTad1File.close()
    outputTad2File.close()

# Termination
command = "rm -rf "+tempLoc
os.system(command)


