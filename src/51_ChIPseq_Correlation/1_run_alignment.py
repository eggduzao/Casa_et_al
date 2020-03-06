
# Import
import os
import sys

# Folder List
indexLoc = "/usr/users/egadegu/Data/Bowtie2Index/hg19.zip"
tempLoc = "/scratch/egadegu/ALG/"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_chip_fastq/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/1_Raw_Bam_Files/"
fastqList = ["SA1_Input", "SA1_Rep1", "SA1_Rep2", "SA2_Input", "SA2_Rep1", "SA2_Rep2"]

# Opening file
inFileName = fl + "1_alg.txt"
inFile = open(inFileName,"w")

# Fasta Loop
for fastqFile in fastqList:

  # Input
  analysisType = "SE"
  minQuality = "20"
  ncores = "4"
  fastqGzFileName = il + fastqFile + ".fastq.gz"
  indexFileName = indexLoc
  tempLocation = tempLoc + fastqFile + "/"
  outputName = fastqFile
  outputLocation = ol

  # Creating files
  inFile.write(" ".join([analysisType, minQuality, ncores, fastqGzFileName, indexFileName, tempLocation, outputName, outputLocation])+"\n")

# Closing file
inFile.close()


