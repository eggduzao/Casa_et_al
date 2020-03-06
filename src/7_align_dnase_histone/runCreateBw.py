
# Import
import os
import sys

# Fastq List
counter = 1
fl = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/7_align_dnase_histone/input_bw/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/"
bamList = ["dnase/bam/HCT116_DNase-seq_UW", "histones/bam/HCT116_ChIP-seq_H3K27ac_USC", "histones/bam/HCT116_ChIP-seq_H3K36me3_USC", "histones/bam/HCT116_ChIP-seq_H3K4me1_USC", "histones/bam/HCT116_ChIP-seq_H3K4me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K4me3_UW", "histones/bam/HCT116_ChIP-seq_H3K9ac_USC", "histones/bam/HCT116_ChIP-seq_H3K9me3_USC", "histones/bam/HCT116_ChIP-seq_H3K9me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K79me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K27me3_BROAD", "histones/bam/HCT116_ChIP-seq_H3K20me1_BROAD", "histones/bam/HCT116_ChIP-seq_H2AFZ_BROAD", "tf/bam/HCT116_ChIP-seq_ATF3_HAIB", "tf/bam/HCT116_ChIP-seq_CBX3_HAIB", "tf/bam/HCT116_ChIP-seq_CEBPB_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_UW", "tf/bam/HCT116_ChIP-seq_EGR1_HAIB", "tf/bam/HCT116_ChIP-seq_ELF1_HAIB", "tf/bam/HCT116_ChIP-seq_EZH2_BROAD", "tf/bam/HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "tf/bam/HCT116_ChIP-seq_FOSL1_HAIB", "tf/bam/HCT116_ChIP-seq_MAX_HAIB", "tf/bam/HCT116_ChIP-seq_POL2RA_USC", "tf/bam/HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "tf/bam/HCT116_ChIP-seq_REST_HAIB", "tf/bam/HCT116_ChIP-seq_SIN3A_HAIB", "tf/bam/HCT116_ChIP-seq_SP1_HAIB", "tf/bam/HCT116_ChIP-seq_TCF7L2_USC", "tf/bam/HCT116_ChIP-seq_TEAD4_HAIB", "tf/bam/HCT116_ChIP-seq_USF1_HAIB", "tf/bam/HCT116_ChIP-seq_YY1_HAIB", "tf/bam/HCT116_ChIP-seq_ZBTB33_HAIB", "tf/bam/HCT116_ChIP-seq_ZNF274_USC", "tf/bam/HCT116_ChIP-seq_SRF_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_BROAD"]

# Fastq Loop
for bamName in bamList:

  # Auxiliary
  name = bamName.split("/")[-1]
  if("tf/bam" in bamName):
    de = "200"
    ue = "0"
    ol = "tf/bw/"
  elif("histones/bam" in bamName):
    de = "200"
    ue = "0"
    ol = "histones/bw/"
  elif("dnase/bam" in bamName):
    de = "5"
    ue = "5"
    ol = "dnase/bw/"

  # Parameters
  downstreamExt = de
  upstreamExt = ue
  genomeSizesFileName = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
  bamFileName = il+bamName+".bam"
  bwFileName = il+ol+name+".bw"

  # Execution
  inputFileName = fl+str(counter)+".txt"
  inputFile = open(inputFileName, "w")
  inputFile.write("\n".join([downstreamExt, upstreamExt, genomeSizesFileName, bamFileName, bwFileName]))
  inputFile.close()
  counter += 1

  # Execution
  #command = "python createBigWig.py "+" ".join([downstreamExt, upstreamExt, genomeSizesFileName, bamFileName, bwFileName])
  #os.system(command)


