
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import itertools
import numpy as np
from functools import reduce

# Input
stagRegionFileName = sys.argv[1]
loopFileName = sys.argv[2]
tempLocation = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def create_anchor_files(inFileName, outFile1Name, outFile2Name):
  
  inFile = open(inFileName, "rU")
  outFile1 = open(outFile1Name, "w")
  outFile2 = open(outFile2Name, "w")
  for line in inFile:
    ll = line.strip().split("\t")
    a1 = int(ll[1]); a2 = int(ll[2])
    outFile1.write("\t".join([ll[0], str(a1), str(a1+25000), ":".join(ll)])+"\n")
    outFile2.write("\t".join([ll[0], str(a2), str(a2+25000), ":".join(ll)])+"\n")
  inFile.close()
  outFile1.close()
  outFile2.close()

def write_anchors_to_file(inFile1Name, inFile2Name, temporary_location, outputFileName):
  
  tempFile1Name = temporary_location + "tempFile1Name.bed"
  inFile1 = open(inFile1Name, "rU")
  inFile2 = open(inFile2Name, "rU")
  tempFile1 = open(tempFile1Name, "w")
  for line in inFile1:
    ll = line.strip().split("\t")
    vec = ll[3].split(":")
    tempFile1.write("\t".join(vec)+"\n")
  for line in inFile2:
    ll = line.strip().split("\t")
    vec = ll[3].split(":")
    tempFile1.write("\t".join(vec)+"\n")
  inFile1.close()
  inFile2.close()
  tempFile1.close()

  tempFile2Name = temporary_location + "tempFile2Name.bed"  
  command = "sort "+tempFile1Name+" | uniq > "+tempFile2Name
  os.system(command)

  command = "sort -k1,1 -k2,2n "+tempFile2Name+" > "+outputFileName
  os.system(command)

def create_intersect_file(stagRegionFileName, loopFileName, tempLocation, outputFileName):

  # Cut files
  cutLoopFileName = tempLocation + "cutLoopFileName.bed"
  cutStagRegionFileName = tempLocation + "cutStagRegionFileName.bed"
  command = "cut -f 1,2,3 "+loopFileName+" > "+cutLoopFileName
  os.system(command)
  command = "cut -f 1,2,3 "+stagRegionFileName+" > "+cutStagRegionFileName
  os.system(command)

  # Creating anchor files
  anchorFile1Name = tempLocation + "anchorFile1Name.bed"
  anchorFile2Name = tempLocation + "anchorFile2Name.bed"
  create_anchor_files(cutLoopFileName, anchorFile1Name, anchorFile2Name)

  # Sorting anchor files and stag file
  sortAnchorFile1Name = tempLocation + "sortAnchorFile1Name.bed"
  sortAnchorFile2Name = tempLocation + "sortAnchorFile2Name.bed"
  sortCutStagRegionFileName = tempLocation + "sortCutStagRegionFileName.bed"
  command = "sort -k1,1 -k2,2n "+anchorFile1Name+" > "+sortAnchorFile1Name
  os.system(command)
  command = "sort -k1,1 -k2,2n "+anchorFile2Name+" > "+sortAnchorFile2Name
  os.system(command)
  command = "sort -k1,1 -k2,2n "+cutStagRegionFileName+" > "+sortCutStagRegionFileName
  os.system(command)

  # Intersecting anchor files with stag file
  intAnchorFile1Name = tempLocation + "intAnchorFile1Name.bed"
  intAnchorFile2Name = tempLocation + "intAnchorFile2Name.bed"
  command = "intersectBed -wa -u -a "+sortAnchorFile1Name+" -b "+sortCutStagRegionFileName+" > "+intAnchorFile1Name
  os.system(command)
  command = "intersectBed -wa -u -a "+sortAnchorFile2Name+" -b "+sortCutStagRegionFileName+" > "+intAnchorFile2Name
  os.system(command)

  # Write anchors to file
  write_anchors_to_file(intAnchorFile1Name, intAnchorFile2Name, tempLocation, outputFileName)

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
create_intersect_file(stagRegionFileName, loopFileName, tempLocation, outputFileName)


