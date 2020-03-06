
# Import
import os
import sys

# Bam List
counter = 1
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/input/"
il = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/1_Raw_Bam_Files/"
tl = "/scratch/egusmao/STAG_MEB/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/2_Merged_Bam_Files/"
bamListList = [["SA1_MCF10A_CHIPSEQ"], ["SA2_MCF10A_CHIPSEQ"], ["SMC1_MCF10A_CHIPSEQ"], ["ZMYM2_MCF10A_CHIPSEQ"], ["Input_MCF10A_CHIPSEQ"], ["INPUT_MCF10A_CHIPSEQ_siSA2"], ["INPUT_MCF10A_CHIPSEQ_siSA1"], ["INPUT_MCF10A_CHIPSEQ_Control"], ["SA1_MCF10A_CHIPSEQ_Control_Rep_1", "SA1_MCF10A_CHIPSEQ_Control_Rep_2"], ["SA2_MCF10A_CHIPSEQ_Control_Rep_1", "SA2_MCF10A_CHIPSEQ_Control_Rep_2"], ["SA1_MCF10A_CHIPSEQ_siSA1_Rep_1", "SA1_MCF10A_CHIPSEQ_siSA1_Rep_2"], ["SA2_MCF10A_CHIPSEQ_siSA1_Rep_1", "SA2_MCF10A_CHIPSEQ_siSA1_Rep_2"], ["SA1_MCF10A_CHIPSEQ_siSA2_Rep_1", "SA1_MCF10A_CHIPSEQ_siSA2_Rep_2"], ["SA2_MCF10A_CHIPSEQ_siSA2_Rep_1", "SA2_MCF10A_CHIPSEQ_siSA2_Rep_2"], ["jc261_MCF10A_CTCF_rep1_CRI01", "jc294_MCF10A_CTCF_rep2_CRI01"], ["jc263_MCF10A_input_CRI01"]]

# Opening file
inFileName = fl + "2_meb.txt"
inFile = open(inFileName,"w")

# Bam Loop
for bamList in bamListList:

  # Parameters
  if("_CRI01" in bamList[0]):
    if(len(bamList) > 1): bamName = "CTCF_MCF10A_CHIPSEQ"
    else: bamName = "INPUT_CTCF_MCF10A_CHIPSEQ"
  elif(len(bamList) == 1): bamName = bamList[0]
  else: bamName = bamList[0].split("_Rep")[0]

  # Input
  bamFileNameList = ",".join([il + e + ".bam" for e in bamList])
  tempBamFileName = tl + str(counter) + "/" + bamName + ".bam"
  outputBamFileName = ol + bamName + ".bam"

  # Creating files
  inFile.write(" ".join([bamFileNameList, tempBamFileName, outputBamFileName])+"\n")
  counter += 1

# Closing file
inFile.close()

