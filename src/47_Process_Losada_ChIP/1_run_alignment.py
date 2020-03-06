
# Import
import os
import sys

# Folder List
indexLoc = "/usr/users/egadegu/Data/Bowtie2Index/hg19.zip"
tempLoc = "/scratch/egadegu/ALG/"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/18_Process_Losada_ChIP/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Data/losada_new/chip/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/1_Raw_Bam_Files/"
fastqList = ["Input_MCF10A_ChIP-seq", "INPUT_MCF10A_ChIP-seq_Control", "INPUT_MCF10A_ChIP-seq_siSA1", "INPUT_MCF10A_ChIP-seq_siSA2", "SA1_MCF10A_ChIP-seq", "SA1_MCF10A_ChIP-seq_Control_Rep_1", "SA1_MCF10A_ChIP-seq_Control_Rep_2", "SA1_MCF10A_ChIP-seq_siSA1_Rep_1", "SA1_MCF10A_ChIP-seq_siSA1_Rep_2", "SA1_MCF10A_ChIP-seq_siSA2_Rep_1", "SA1_MCF10A_ChIP-seq_siSA2_Rep_2", "SA2_MCF10A_ChIP-seq", "SA2_MCF10A_ChIP-seq_Control_Rep_1", "SA2_MCF10A_ChIP-seq_Control_Rep_2", "SA2_MCF10A_ChIP-seq_siSA1_Rep_1", "SA2_MCF10A_ChIP-seq_siSA1_Rep_2", "SA2_MCF10A_ChIP-seq_siSA2_Rep_1", "SA2_MCF10A_ChIP-seq_siSA2_Rep_2", "SMC1_MCF10A_ChIP-seq", "ZMYM2_MCF10A_ChIP-seq"]

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


