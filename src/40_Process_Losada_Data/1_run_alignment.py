
# Import
import os
import sys

# Folder List
counter = 1
indexLoc = "/projects/ag-papan/genomes/BowtieIndexes/hg19.zip"
tempLoc = "/scratch/egusmao/STAG_ALG/"
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/input/"
il = "/projects/ag-papan/eduardo/Wendt_Stag/Data/losada/chip/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/1_Raw_Bam_Files/"
fastqList = ["SA1_MCF10A_CHIPSEQ", "SA2_MCF10A_CHIPSEQ", "SMC1_MCF10A_CHIPSEQ", "ZMYM2_MCF10A_CHIPSEQ", "Input_MCF10A_CHIPSEQ", "SA1_MCF10A_CHIPSEQ_Control_Rep_1", "SA1_MCF10A_CHIPSEQ_Control_Rep_2", "SA2_MCF10A_CHIPSEQ_Control_Rep_1", "SA2_MCF10A_CHIPSEQ_Control_Rep_2", "SA1_MCF10A_CHIPSEQ_siSA1_Rep_1", "SA1_MCF10A_CHIPSEQ_siSA1_Rep_2", "SA2_MCF10A_CHIPSEQ_siSA1_Rep_1", "SA2_MCF10A_CHIPSEQ_siSA1_Rep_2", "SA1_MCF10A_CHIPSEQ_siSA2_Rep_1", "SA1_MCF10A_CHIPSEQ_siSA2_Rep_2", "SA2_MCF10A_CHIPSEQ_siSA2_Rep_1", "SA2_MCF10A_CHIPSEQ_siSA2_Rep_2", "INPUT_MCF10A_CHIPSEQ_siSA2", "INPUT_MCF10A_CHIPSEQ_siSA1", "INPUT_MCF10A_CHIPSEQ_Control", "jc261_MCF10A_CTCF_rep1_CRI01", "jc294_MCF10A_CTCF_rep2_CRI01", "jc263_MCF10A_input_CRI01"]

# Opening file
inFileName = fl + "1_alg.txt"
inFile = open(inFileName,"w")

# Fasta Loop
for fastqFile in fastqList:

  # Input
  analysisType = "SE"
  minQuality = "20"
  ncores = "1"
  fastqGzFileName = il + fastqFile + ".fastq.gz"
  indexFileName = indexLoc
  tempLocation = tempLoc + str(counter) + "/"
  outputName = fastqFile
  outputLocation = ol

  # Creating files
  inFile.write(" ".join([analysisType, minQuality, ncores, fastqGzFileName, indexFileName, tempLocation, outputName, outputLocation])+"\n")
  counter += 1

# Closing file
inFile.close()


