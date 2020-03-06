
#Import
import os
import sys

# Macs list
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/"
macsNameList = ["dnase/macs/HCT116_DNase-seq_UW", "histones/macs/HCT116_ChIP-seq_H3K27ac_USC", "histones/macs/HCT116_ChIP-seq_H3K36me3_USC", "histones/macs/HCT116_ChIP-seq_H3K4me1_USC", "histones/macs/HCT116_ChIP-seq_H3K4me2_BROAD", "histones/macs/HCT116_ChIP-seq_H3K4me3_UW", "histones/macs/HCT116_ChIP-seq_H3K9ac_USC", "histones/macs/HCT116_ChIP-seq_H3K9me3_USC", "histones/macs/HCT116_ChIP-seq_H3K9me2_BROAD", "histones/macs/HCT116_ChIP-seq_H3K79me2_BROAD", "histones/macs/HCT116_ChIP-seq_H3K27me3_BROAD", "histones/macs/HCT116_ChIP-seq_H3K20me1_BROAD", "histones/macs/HCT116_ChIP-seq_H2AFZ_BROAD", "tf/macs/HCT116_ChIP-seq_ATF3_HAIB", "tf/macs/HCT116_ChIP-seq_CBX3_HAIB", "tf/macs/HCT116_ChIP-seq_CEBPB_HAIB", "tf/macs/HCT116_ChIP-seq_CTCF_HAIB", "tf/macs/HCT116_ChIP-seq_CTCF_UW", "tf/macs/HCT116_ChIP-seq_EGR1_HAIB", "tf/macs/HCT116_ChIP-seq_ELF1_HAIB", "tf/macs/HCT116_ChIP-seq_EZH2_BROAD", "tf/macs/HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "tf/macs/HCT116_ChIP-seq_FOSL1_HAIB", "tf/macs/HCT116_ChIP-seq_MAX_HAIB", "tf/macs/HCT116_ChIP-seq_POL2RA_USC", "tf/macs/HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "tf/macs/HCT116_ChIP-seq_REST_HAIB", "tf/macs/HCT116_ChIP-seq_SIN3A_HAIB", "tf/macs/HCT116_ChIP-seq_SP1_HAIB", "tf/macs/HCT116_ChIP-seq_TCF7L2_USC", "tf/macs/HCT116_ChIP-seq_TEAD4_HAIB", "tf/macs/HCT116_ChIP-seq_USF1_HAIB", "tf/macs/HCT116_ChIP-seq_YY1_HAIB", "tf/macs/HCT116_ChIP-seq_ZBTB33_HAIB", "tf/macs/HCT116_ChIP-seq_ZNF274_USC", "tf/macs/HCT116_ChIP-seq_SRF_HAIB", "tf/macs/HCT116_ChIP-seq_CTCF_BROAD", "tf/macs/HCT116_ChIP-seq_JUND_HAIB", "tf/macs/HCT116_ChIP-seq_RAD21_HAIB", "tf/macs/HCT116_ChIP-seq_ZFX_SYDH"]

macsNameList = ["tf/macs/HCT116_ChIP-seq_JUND_HAIB", "tf/macs/HCT116_ChIP-seq_RAD21_HAIB", "tf/macs/HCT116_ChIP-seq_ZFX_SYDH"]

# Macs Loop
for macsName in macsNameList:

  # Input
  peakFileName = il+macsName+"_peaks.narrowPeak"
  summitFileName = il+macsName+"_summits.bed"

  # Execution
  print macsName
  command = "python fixMacsOutput.py "+" ".join([peakFileName, summitFileName])
  os.system(command)


