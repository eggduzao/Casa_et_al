
# Import
import os
import sys

# Hic List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/19_Process_Losada_Hic/input/"
juicerCommand = "/usr/users/egadegu/Juicer/scripts/common/mega.sh"
genomeID = "hg19"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/19_Process_Losada_Hic/1_Juicer/"
tl = "/scratch/egadegu/RMH/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/19_Process_Losada_Hic/2_Merged_Hic/"
hicFolderList = [["LCONT_rep1", "LCONT_rep2"], ["SISA1_rep1", "SISA1_rep2"], ["SISA2_rep1", "SISA2_rep2"]]

# Opening file
inFileName = fl + "2_rmh.txt"
inFile = open(inFileName,"w")

# Hic Loop
for hicFolderVec in hicFolderList:

  # Parameters
  finalName = hicFolderVec[0].split("_")[0]

  # Input
  juicerCommandFile = juicerCommand
  genomeId = genomeID
  restrictionEnzyme = "none"
  excludeFragments = "1"
  hicWorkFolderList = ",".join([il + e + "/" for e in hicFolderVec])
  tempLocation = tl + finalName + "/"
  outputLocation = ol + finalName + "/"

  # Creating files
  inFile.write(" ".join([juicerCommandFile, genomeId, restrictionEnzyme, excludeFragments, hicWorkFolderList, tempLocation, outputLocation])+"\n")

inFile.close()


