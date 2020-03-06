
# Import
import os
import sys

# Folder List
indexLoc = "/usr/users/egadegu/Data/Bowtie2Index/hg19.zip"
tempLoc = "/scratch/egadegu/ALG/"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Data/AKI69_Wendt_Smc3_ChIPseq/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/1_Raw_Bam_Files/"
fastqList = ["Smc3_Input", "Smc3_Rep1", "Smc3_Rep2"]

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


