
# Import
import os
import sys
from pysam import Samfile

# Type List
regionsFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/1_genomic_regions/regions.bam"
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/2_stag_regions_overlap/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/3_barplot/"
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
  # Without SE / Without active status
  #####################################

  # Input
  legendRows = "5"
  inputStag1FileName = tempLocation+"inputStag1FileName.txt"
  inputStag2FileName = tempLocation+"inputStag2FileName.txt"
  outputBarFileName = ol+"regions_woSE_woST.pdf"
  outputTotalBarFileName = ol+"totalSignal_woSE_woST.pdf"
  outputViolinFileName = ol+"violinSignal_woSE_woST.pdf"
  outputTotalCtcfBarFileName = ol+"totalCtcfSignal_woSE_woST.pdf"
  outputCtcfViolinFileName = ol+"violinCtcfSignal_woSE_woST.pdf"

  # Iterating on STAG1
  stag1File = open(stag1FileName, "rU")
  inputStag1File = open(inputStag1FileName, "w")
  stag1File.readline()
  for line in stag1File:
    ll = line.strip().split("\t")
    if(ll[7].split(":")[0] == "SUPERENHANCER"): continue
    name = ll[7].split(":")[0]
    score = ll[3]
    ctcfScore = ll[4]
    inputStag1File.write("\t".join([name, score, ctcfScore])+"\n")
  stag1File.close()
  inputStag1File.close()

  # Iterating on STAG2
  stag2File = open(stag2FileName, "rU")
  inputStag2File = open(inputStag2FileName, "w")
  stag2File.readline()
  for line in stag2File:
    ll = line.strip().split("\t")
    if(ll[7].split(":")[0] == "SUPERENHANCER"): continue
    name = ll[7].split(":")[0]
    score = ll[3]
    ctcfScore = ll[4]
    inputStag2File.write("\t".join([name, score, ctcfScore])+"\n")
  stag2File.close()
  inputStag2File.close()

  # Execution
  command = "R CMD BATCH '--args '"+legendRows+"' '"+inputStag1FileName+"' '"+inputStag2FileName+"' '"+outputBarFileName+"' '"+outputTotalBarFileName+"' '"+outputViolinFileName+"' '"+outputTotalCtcfBarFileName+"' '"+outputCtcfViolinFileName+" 3_barplot.R 3_barplot.Rout"
  os.system(command)  

  #####################################
  # Without SE / With active status
  #####################################

  # Input
  legendRows = "8"
  inputStag1FileName = tempLocation+"inputStag1FileName.txt"
  inputStag2FileName = tempLocation+"inputStag2FileName.txt"
  outputBarFileName = ol+"regions_woSE_withST.pdf"
  outputTotalBarFileName = ol+"totalSignal_woSE_withST.pdf"
  outputViolinFileName = ol+"violinSignal_woSE_withST.pdf"
  outputTotalCtcfBarFileName = ol+"totalCtcfSignal_woSE_withST.pdf"
  outputCtcfViolinFileName = ol+"violinCtcfSignal_woSE_withST.pdf"

  # Iterating on STAG1
  stag1File = open(stag1FileName, "rU")
  inputStag1File = open(inputStag1FileName, "w")
  stag1File.readline()
  for line in stag1File:
    ll = line.strip().split("\t")
    sp = ll[7].split(":")
    if(sp[0] == "SUPERENHANCER"): continue
    if(sp[2] == "."): name = sp[0]
    else: name = sp[2]+"_"+sp[0]
    score = ll[3]
    ctcfScore = ll[4]
    inputStag1File.write("\t".join([name, score, ctcfScore])+"\n")
  stag1File.close()
  inputStag1File.close()

  # Iterating on STAG2
  stag2File = open(stag2FileName, "rU")
  inputStag2File = open(inputStag2FileName, "w")
  stag2File.readline()
  for line in stag2File:
    ll = line.strip().split("\t")
    sp = ll[7].split(":")
    if(sp[0] == "SUPERENHANCER"): continue
    if(sp[2] == "."): name = sp[0]
    else: name = sp[2]+"_"+sp[0]
    score = ll[3]
    ctcfScore = ll[4]
    inputStag2File.write("\t".join([name, score, ctcfScore])+"\n")
  stag2File.close()
  inputStag2File.close()

  # Execution
  command = "R CMD BATCH '--args '"+legendRows+"' '"+inputStag1FileName+"' '"+inputStag2FileName+"' '"+outputBarFileName+"' '"+outputTotalBarFileName+"' '"+outputViolinFileName+"' '"+outputTotalCtcfBarFileName+"' '"+outputCtcfViolinFileName+" 3_barplot.R 3_barplot.Rout"
  os.system(command)

  #####################################
  # With SE / Without active status
  #####################################

  # Input
  legendRows = "6"
  inputStag1FileName = tempLocation+"inputStag1FileName.txt"
  inputStag2FileName = tempLocation+"inputStag2FileName.txt"
  outputBarFileName = ol+"regions_withSE_woST.pdf"
  outputTotalBarFileName = ol+"totalSignal_withSE_woST.pdf"
  outputViolinFileName = ol+"violinSignal_withSE_woST.pdf"
  outputTotalCtcfBarFileName = ol+"totalCtcfSignal_withSE_woST.pdf"
  outputCtcfViolinFileName = ol+"violinCtcfSignal_withSE_woST.pdf"

  # Iterating on STAG1
  stag1File = open(stag1FileName, "rU")
  inputTempStag1File = tempLocation+"inputStag1FileName.txt"
  inputStag1File = open(inputStag1FileName, "w")
  stag1File.readline()
  for line in stag1File:
    ll = line.strip().split("\t")
    name = ll[7].split(":")[0]
    if(name == "ENHANCER"):
      regionsFetch = regionsFile.fetch(ll[0], int(ll[5]), int(ll[6]))
      seFlag = False
      for read in regionsFetch:
        if(read.qname.split(":")[0] == "SUPERENHANCER"):
          seFlag = True
          break
      if(seFlag): continue
    score = ll[3]
    ctcfScore = ll[4]
    inputStag1File.write("\t".join([name, score, ctcfScore])+"\n")
  stag1File.close()
  inputStag1File.close()

  # Iterating on STAG2
  stag2File = open(stag2FileName, "rU")
  inputStag2File = open(inputStag2FileName, "w")
  stag2File.readline()
  for line in stag2File:
    ll = line.strip().split("\t")
    name = ll[7].split(":")[0]
    if(name == "ENHANCER"):
      regionsFetch = regionsFile.fetch(ll[0], int(ll[5]), int(ll[6]))
      seFlag = False
      for read in regionsFetch:
        if(read.qname.split(":")[0] == "SUPERENHANCER"):
          seFlag = True
          break
      if(seFlag): continue
    score = ll[3]
    ctcfScore = ll[4]
    inputStag2File.write("\t".join([name, score, ctcfScore])+"\n")
  stag2File.close()
  inputStag2File.close()

  # Execution
  command = "R CMD BATCH '--args '"+legendRows+"' '"+inputStag1FileName+"' '"+inputStag2FileName+"' '"+outputBarFileName+"' '"+outputTotalBarFileName+"' '"+outputViolinFileName+"' '"+outputTotalCtcfBarFileName+"' '"+outputCtcfViolinFileName+" 3_barplot.R 3_barplot.Rout"
  os.system(command) 

  #####################################
  # With SE / With active status
  #####################################

  # Input
  legendRows = "10"
  inputStag1FileName = tempLocation+"inputStag1FileName.txt"
  inputStag2FileName = tempLocation+"inputStag2FileName.txt"
  outputBarFileName = ol+"regions_withSE_withST.pdf"
  outputTotalBarFileName = ol+"totalSignal_withSE_withST.pdf"
  outputViolinFileName = ol+"violinSignal_withSE_withST.pdf"
  outputTotalCtcfBarFileName = ol+"totalCtcfSignal_withSE_withST.pdf"
  outputCtcfViolinFileName = ol+"violinCtcfSignal_withSE_withST.pdf"

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
    inputStag1File.write("\t".join([name, score, ctcfScore])+"\n")
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
    inputStag2File.write("\t".join([name, score, ctcfScore])+"\n")
  stag2File.close()
  inputStag2File.close()

  # Execution
  command = "R CMD BATCH '--args '"+legendRows+"' '"+inputStag1FileName+"' '"+inputStag2FileName+"' '"+outputBarFileName+"' '"+outputTotalBarFileName+"' '"+outputViolinFileName+"' '"+outputTotalCtcfBarFileName+"' '"+outputCtcfViolinFileName+" 3_barplot.R 3_barplot.Rout"
  os.system(command)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)

