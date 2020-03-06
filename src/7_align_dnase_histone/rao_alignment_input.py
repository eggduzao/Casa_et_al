
# Import
import os
import sys

# Fastq List
counter = 1
fl = "./rao_alignment/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/rao/fasta/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/rao/bam/"
fastqList = ["CHIP001_RAD21_untreated_SE_1_1", "CHIP001_RAD21_untreated_SE_2_1", "CHIP001_RAD21_untreated_SE_3_1", "CHIP001_RAD21_untreated_SE_4_1", "CHIP002_RAD21_treated_SE_1_1", "CHIP002_RAD21_treated_SE_2_1", "CHIP002_RAD21_treated_SE_3_1", "CHIP002_RAD21_treated_SE_4_1", ["CHIP003_SMC1_untreated_PE_1_1", "CHIP003_SMC1_untreated_PE_1_2"], ["CHIP003_SMC1_untreated_PE_2_1", "CHIP003_SMC1_untreated_PE_2_2"], ["CHIP003_SMC1_untreated_PE_3_1", "CHIP003_SMC1_untreated_PE_3_2"], ["CHIP003_SMC1_untreated_PE_4_1", "CHIP003_SMC1_untreated_PE_4_2"], ["CHIP004_SMC1_treated_PE_1_1", "CHIP004_SMC1_treated_PE_1_2"], ["CHIP004_SMC1_treated_PE_2_1", "CHIP004_SMC1_treated_PE_2_2"], ["CHIP004_SMC1_treated_PE_3_1", "CHIP004_SMC1_treated_PE_3_2"], ["CHIP004_SMC1_treated_PE_4_1", "CHIP004_SMC1_treated_PE_4_2"], ["CHIP005_CTCF_untreated_PE_1_1", "CHIP005_CTCF_untreated_PE_1_2"], ["CHIP005_CTCF_untreated_PE_2_1", "CHIP005_CTCF_untreated_PE_2_2"], ["CHIP005_CTCF_untreated_PE_3_1", "CHIP005_CTCF_untreated_PE_3_2"], ["CHIP005_CTCF_untreated_PE_4_1", "CHIP005_CTCF_untreated_PE_4_2"], ["CHIP006_CTCF_treated_PE_1_1", "CHIP006_CTCF_treated_PE_1_2"], ["CHIP006_CTCF_treated_PE_2_1", "CHIP006_CTCF_treated_PE_2_2"], ["CHIP006_CTCF_treated_PE_3_1", "CHIP006_CTCF_treated_PE_3_2"], ["CHIP006_CTCF_treated_PE_4_1", "CHIP006_CTCF_treated_PE_4_2"], "CHIP007_CTCF_untreated_SE_1_1", "CHIP007_CTCF_untreated_SE_2_1", "CHIP007_CTCF_untreated_SE_3_1", "CHIP007_CTCF_untreated_SE_4_1", "CHIP008_CTCF_treated_SE_1_1", "CHIP008_CTCF_treated_SE_2_1", "CHIP009_H3K27Ac_untreated_SE_1_1", "CHIP010_H3K27Ac_treated_SE_1_1", "CHIP011_H3K4me3_untreated_SE_1_1", "CHIP012_H3K4me3_treated_SE_1_1", "CHIP013_H3K4me1_untreated_SE_1_1", "CHIP014_H3K4me1_treated_SE_1_1", "CHIP015_H3K36me3_untreated_SE_1_1", "CHIP016_H3K36me3_treated_SE_1_1", "CHIP017_H3K27me3_untreated_SE_1_1", "CHIP018_H3K27me3_treated_SE_1_1", "CHIP019_H3K9me3_untreated_SE_1_1", "CHIP020_H3K9me3_treated_SE_1_1", "CHIP021_H4K16Ac_untreated_SE_1_1", "CHIP022_H4K16Ac_treated_SE_1_1", "CHIP023_H3K79me2_untreated_SE_1_1", "CHIP024_H3K79me2_treated_SE_1_1", "CHIP025_H4K20me3_untreated_SE_1_1", "CHIP026_H4K20me3_treated_SE_1_1", "CHIP027_H2AZ_untreated_SE_1_1", "CHIP028_H2AZ_treated_SE_1_1", "CHIP029_NIPBL_untreated_SE_1_1", "CHIP029_NIPBL_untreated_SE_2_1", "CHIP029_NIPBL_untreated_SE_3_1", "CHIP029_NIPBL_untreated_SE_4_1", "CHIP030_NIPBL_treated_SE_1_1", "CHIP031_input_untreated_SE_1_1", "CHIP031_input_untreated_SE_2_1", "CHIP031_input_untreated_SE_3_1", "CHIP031_input_untreated_SE_4_1", "CHIP032_input_treated_SE_1_1", "CHIP032_input_treated_SE_2_1", "CHIP032_input_treated_SE_3_1", "CHIP032_input_treated_SE_4_1", ["CHIP033_input_untreated_PE_1_1", "CHIP033_input_untreated_PE_1_2"], ["CHIP033_input_untreated_PE_2_1", "CHIP033_input_untreated_PE_2_2"], ["CHIP033_input_untreated_PE_3_1", "CHIP033_input_untreated_PE_3_2"], ["CHIP033_input_untreated_PE_4_1", "CHIP033_input_untreated_PE_4_2"], ["CHIP034_input_treated_PE_1_1", "CHIP034_input_treated_PE_1_2"], ["CHIP034_input_treated_PE_2_1", "CHIP034_input_treated_PE_2_2"], ["CHIP034_input_treated_PE_3_1", "CHIP034_input_treated_PE_3_2"], ["CHIP034_input_treated_PE_4_1", "CHIP034_input_treated_PE_4_2"], "CHIP035_input_untreated_SE_1_1", "CHIP036_input_treated_SE_1_1", "CHIP037_input_untreated_SE_1_1", "CHIP038_input_treated_SE_1_1", "CHIP039_input_treated_SE_1_1"]

# Fastq Loop
for fastq in fastqList:

  # SE or PE
  if(isinstance(fastq, basestring)):

    # Parameters
    alignType = "SE"
    minQuality = "20"
    ncores = "1"
    fastqFileName = il + fastq + ".fastq"
    indexFileName = "/projects/ag-papan/genomes/BowtieIndexes/hg19.zip"
    tempLocation = "/scratch/eduardo/stag_align_rao_se/"
    outputLocation = ol

    # Creating files
    inFileName = fl+str(counter)+".txt"
    inFile = open(inFileName,"w")
    inFile.write("\n".join([alignType, minQuality, ncores, fastqFileName, indexFileName, tempLocation, outputLocation]))
    inFile.close()
    counter += 1

  elif(isinstance(fastq, list)):

    # Parameters
    alignType = "PE"
    minQuality = "20"
    ncores = "1"
    fastqFileName = il + fastq[0] + ".fastq," + il + fastq[1] + ".fastq"
    indexFileName = "/projects/ag-papan/genomes/BowtieIndexes/hg19.zip"
    tempLocation = "/scratch/eduardo/stag_align_rao_pe/"
    outputLocation = ol

    # Creating files
    inFileName = fl+str(counter)+".txt"
    inFile = open(inFileName,"w")
    inFile.write("\n".join([alignType, minQuality, ncores, fastqFileName, indexFileName, tempLocation, outputLocation]))
    inFile.close()
    counter += 1


