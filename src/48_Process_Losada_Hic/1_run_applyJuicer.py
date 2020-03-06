
# Import
import os
import sys

# Folder List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/19_Process_Losada_Hic/input/"
juicerCommand = "/usr/users/egadegu/Juicer/scripts/juicer.sh"
genomeID = "hg19"
chromSizesFile = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19"
rl = "/usr/users/egadegu/Juicer/restriction_sites/hg19_"
genomeFile = "/usr/users/egadegu/Juicer/references/Homo_sapiens_assembly19.fasta"
juicerLoc = "/usr/users/egadegu/Juicer/"
wl = "/scratch/egadegu/WENDT_LOSADA_JUICER/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Data/losada_new/hic/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/19_Process_Losada_Hic/1_Juicer/"
fastqList = [["LCONT_rep1_R1", "LCONT_rep1_R2"], ["LCONT_rep2_R1", "LCONT_rep2_R2"], ["SISA1_rep1_R1", "SISA1_rep1_R2"], ["SISA1_rep2_R1", "SISA1_rep2_R2"], ["SISA2_rep1_R1", "SISA2_rep1_R2"], ["SISA2_rep2_R1", "SISA2_rep2_R2"]]

# Opening file
inFileName = fl + "1_jui.txt"
inFile = open(inFileName,"w")

# Fasta Loop
for fastqPair in fastqList:

  # Parameters
  folderName = "_".join(fastqPair[0].split("_")[:-1])
  expName = "_".join(fastqPair[0].split("_")[:-1])

  # Input
  juicerCommandFile = juicerCommand
  genomeId = genomeID
  restrictionEnzyme = "none"
  expDescription = expName
  commandStage = "."
  useShort = "0"
  chromSizesLocation = chromSizesFile
  restrictionSiteFileName = rl + "none" + ".txt"
  assemblyGenomeFileName = genomeFile
  juicerFolder = juicerLoc
  numberOfThreads = "20"
  fastq1FileName = il + fastqPair[0] + ".fastq.gz"
  fastq2FileName = il + fastqPair[1] + ".fastq.gz"
  workingFolder = wl + expName + "/"
  outputLocation = ol + folderName + "/" + expName + "/"

  # Creating files
  inFile.write(" ".join([juicerCommandFile, genomeId, restrictionEnzyme, expDescription, commandStage, useShort, chromSizesLocation, restrictionSiteFileName, assemblyGenomeFileName, juicerFolder, numberOfThreads, fastq1FileName, fastq2FileName, workingFolder, outputLocation])+"\n")

# Closing file
inFile.close()




