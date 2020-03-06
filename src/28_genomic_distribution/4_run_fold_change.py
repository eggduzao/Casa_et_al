
# Import
import os
import sys
from pysam import Samfile

# Type List
regionsFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/1_genomic_regions/regions.bam"
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/2_stag_regions_overlap/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/4_fold_change/"
typeList = ["_regions.txt"]

# Initialization
tempLocation = ol+"TEMP/"
command = "mkdir -p "+tempLocation
os.system(command)
regionsFile = Samfile(regionsFileName, "rb")

# Type Loop
for typ in typeList:

  # Input files
  stag1FileName = il+"STAG1"+typ
  stag2FileName = il+"STAG2"+typ

  #####################################
  # With SE / With active status
  #####################################

  # Input
  inputStag1FileName = tempLocation+"inputStag1FileName.txt"
  inputStag2FileName = tempLocation+"inputStag2FileName.txt"
  outputRegionFileName = ol+"region_fold_change.pdf"
  outputSignalFileName = ol+"signal_fold_change.pdf"

  # Iterating on STAG1
  stag1File = open(stag1FileName, "rU")
  inputStag1File = open(inputStag1FileName, "w")
  stag1File.readline()
  for line in stag1File:
    ll = line.strip().split("\t")
    sp = ll[7].split(":")
    if(sp[0] == "ENHANCER"):
      regionsFetch = regionsFile.fetch(ll[0], int(ll[5]), int(ll[6]))
      seFlag = False
      for read in regionsFetch:
        if(read.qname.split(":")[0] == "SUPERENHANCER"):
          seFlag = True
          break
      if(seFlag): continue
    if(sp[2] == "."): name = sp[0]
    else: name = sp[2]+"_"+sp[0]
    score = ll[3]
    ctcfScore = ll[4]
    inputStag1File.write("\t".join([name, score])+"\n")
  stag1File.close()
  inputStag1File.close()

  # Iterating on STAG2
  stag2File = open(stag2FileName, "rU")
  inputStag2File = open(inputStag2FileName, "w")
  stag2File.readline()
  for line in stag2File:
    ll = line.strip().split("\t")
    sp = ll[7].split(":")
    if(sp[0] == "ENHANCER"):
      regionsFetch = regionsFile.fetch(ll[0], int(ll[5]), int(ll[6]))
      seFlag = False
      for read in regionsFetch:
        if(read.qname.split(":")[0] == "SUPERENHANCER"):
          seFlag = True
          break
      if(seFlag): continue
    if(sp[2] == "."): name = sp[0]
    else: name = sp[2]+"_"+sp[0]
    score = ll[3]
    ctcfScore = ll[4]
    inputStag2File.write("\t".join([name, score])+"\n")
  stag2File.close()
  inputStag2File.close()

  # Execution
  command = "R CMD BATCH '--args '"+inputStag1FileName+"' '"+inputStag2FileName+"' '"+outputRegionFileName+"' '"+outputSignalFileName+" 4_foldChange.R 4_foldChange.Rout"
  os.system(command)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)

