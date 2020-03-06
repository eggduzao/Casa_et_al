
# Import
import os
import sys

# Fastq List
counter = 1
fl = "./input_macs/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/"
treatList = ["dnase/bam/HCT116_DNase-seq_UW", "histones/bam/HCT116_ChIP-seq_H3K27ac_USC", "histones/bam/HCT116_ChIP-seq_H3K36me3_USC", "histones/bam/HCT116_ChIP-seq_H3K4me1_USC", "histones/bam/HCT116_ChIP-seq_H3K4me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K4me3_UW", "histones/bam/HCT116_ChIP-seq_H3K9ac_USC", "histones/bam/HCT116_ChIP-seq_H3K9me3_USC", "histones/bam/HCT116_ChIP-seq_H3K9me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K79me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K27me3_BROAD", "histones/bam/HCT116_ChIP-seq_H3K20me1_BROAD", "histones/bam/HCT116_ChIP-seq_H2AFZ_BROAD", "tf/bam/HCT116_ChIP-seq_ATF3_HAIB", "tf/bam/HCT116_ChIP-seq_CBX3_HAIB", "tf/bam/HCT116_ChIP-seq_CEBPB_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_UW", "tf/bam/HCT116_ChIP-seq_EGR1_HAIB", "tf/bam/HCT116_ChIP-seq_ELF1_HAIB", "tf/bam/HCT116_ChIP-seq_EZH2_BROAD", "tf/bam/HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "tf/bam/HCT116_ChIP-seq_FOSL1_HAIB", "tf/bam/HCT116_ChIP-seq_MAX_HAIB", "tf/bam/HCT116_ChIP-seq_POL2RA_USC", "tf/bam/HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "tf/bam/HCT116_ChIP-seq_REST_HAIB", "tf/bam/HCT116_ChIP-seq_SIN3A_HAIB", "tf/bam/HCT116_ChIP-seq_SP1_HAIB", "tf/bam/HCT116_ChIP-seq_TCF7L2_USC", "tf/bam/HCT116_ChIP-seq_TEAD4_HAIB", "tf/bam/HCT116_ChIP-seq_USF1_HAIB", "tf/bam/HCT116_ChIP-seq_YY1_HAIB", "tf/bam/HCT116_ChIP-seq_ZBTB33_HAIB", "tf/bam/HCT116_ChIP-seq_ZNF274_USC", "tf/bam/hg19_HCT116_ChIP-seq_SRF_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_BROAD", "tf/bam/HCT116_ChIP-seq_JUND_HAIB", "tf/bam/HCT116_ChIP-seq_RAD21_HAIB", "tf/bam/HCT116_ChIP-seq_ZFX_SYDH"]

# Fastq Loop
for treat in treatList:

  # Auxiliary
  if("_ZFX_" in treat):
    dtype = "PE"
    ol = "tf/macs/"
  elif("tf/bam" in treat):
    dtype = "ChIP-seq"
    ol = "tf/macs/"
  elif("histones/bam" in treat):
    dtype = "ChIP-seq"
    ol = "histones/macs/"
  elif("dnase/bam" in treat):
    dtype = "DNase-seq"
    ol = "dnase/macs/"

  if("dnase/bam" in treat): contr = "."
  elif("_BROAD" in treat): contr = "control/bam/HCT116_ChIP-seq_Control_BROAD"
  elif("_HAIB" in treat): contr = "control/bam/HCT116_ChIP-seq_Control_HAIB"
  elif("_USC" in treat): contr = "control/bam/HCT116_ChIP-seq_Control_USC"
  elif("_UW" in treat): contr = "control/bam/HCT116_ChIP-seq_Control_UW"
  else: contr = "."

  # Parameters
  dataType = dtype
  treatmentFileName = il+treat+".bam"
  controlFileName = il+contr+".bam"
  tempLocation = "/scratch/eduardo/stag_macs/" + treat.split("/")[-1] + "/"
  outputLocation = il+ol

  # Creating files
  inFileName = fl+str(counter)+".txt"
  inFile = open(inFileName,"w")
  inFile.write("\n".join([dataType, treatmentFileName, controlFileName, tempLocation, outputLocation]))
  inFile.close()
  counter += 1


