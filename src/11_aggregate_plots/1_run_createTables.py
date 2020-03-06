
# Import
import os
import sys

# Input List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/"
ol = "/home/egg/Projects/Papantonis_Stag/Results/11_aggregate_plots/tables/"
bamList = ["dnase/bam/HCT116_DNase-seq_UW", "histones/bam/HCT116_ChIP-seq_H3K27ac_USC", "histones/bam/HCT116_ChIP-seq_H3K36me3_USC", "histones/bam/HCT116_ChIP-seq_H3K4me1_USC", "histones/bam/HCT116_ChIP-seq_H3K4me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K4me3_UW", "histones/bam/HCT116_ChIP-seq_H3K9ac_USC", "histones/bam/HCT116_ChIP-seq_H3K9me3_USC", "histones/bam/HCT116_ChIP-seq_H3K9me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K79me2_BROAD", "histones/bam/HCT116_ChIP-seq_H3K27me3_BROAD", "histones/bam/HCT116_ChIP-seq_H3K20me1_BROAD", "histones/bam/HCT116_ChIP-seq_H2AFZ_BROAD", "tf/bam/HCT116_ChIP-seq_ATF3_HAIB", "tf/bam/HCT116_ChIP-seq_CBX3_HAIB", "tf/bam/HCT116_ChIP-seq_CEBPB_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_UW", "tf/bam/HCT116_ChIP-seq_EGR1_HAIB", "tf/bam/HCT116_ChIP-seq_ELF1_HAIB", "tf/bam/HCT116_ChIP-seq_EZH2_BROAD", "tf/bam/HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "tf/bam/HCT116_ChIP-seq_FOSL1_HAIB", "tf/bam/HCT116_ChIP-seq_MAX_HAIB", "tf/bam/HCT116_ChIP-seq_POL2RA_USC", "tf/bam/HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "tf/bam/HCT116_ChIP-seq_REST_HAIB", "tf/bam/HCT116_ChIP-seq_SIN3A_HAIB", "tf/bam/HCT116_ChIP-seq_SP1_HAIB", "tf/bam/HCT116_ChIP-seq_TCF7L2_USC", "tf/bam/HCT116_ChIP-seq_TEAD4_HAIB", "tf/bam/HCT116_ChIP-seq_USF1_HAIB", "tf/bam/HCT116_ChIP-seq_YY1_HAIB", "tf/bam/HCT116_ChIP-seq_ZBTB33_HAIB", "tf/bam/HCT116_ChIP-seq_ZNF274_USC", "tf/bam/hg19_HCT116_ChIP-seq_SRF_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_BROAD"]

bamList = ["tf/bam/HCT116_ChIP-seq_ATF3_HAIB", "tf/bam/HCT116_ChIP-seq_CBX3_HAIB", "tf/bam/HCT116_ChIP-seq_CEBPB_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_UW", "tf/bam/HCT116_ChIP-seq_EGR1_HAIB", "tf/bam/HCT116_ChIP-seq_ELF1_HAIB", "tf/bam/HCT116_ChIP-seq_EZH2_BROAD", "tf/bam/HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "tf/bam/HCT116_ChIP-seq_FOSL1_HAIB", "tf/bam/HCT116_ChIP-seq_MAX_HAIB", "tf/bam/HCT116_ChIP-seq_POL2RA_USC", "tf/bam/HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "tf/bam/HCT116_ChIP-seq_REST_HAIB", "tf/bam/HCT116_ChIP-seq_SIN3A_HAIB", "tf/bam/HCT116_ChIP-seq_SP1_HAIB", "tf/bam/HCT116_ChIP-seq_TCF7L2_USC", "tf/bam/HCT116_ChIP-seq_TEAD4_HAIB", "tf/bam/HCT116_ChIP-seq_USF1_HAIB", "tf/bam/HCT116_ChIP-seq_YY1_HAIB", "tf/bam/HCT116_ChIP-seq_ZBTB33_HAIB", "tf/bam/HCT116_ChIP-seq_ZNF274_USC", "tf/bam/hg19_HCT116_ChIP-seq_SRF_HAIB", "tf/bam/HCT116_ChIP-seq_CTCF_BROAD"]

# Input Loop
for bamName in bamList:

  # Bed List
  #bl1 = "/home/egg/Projects/Papantonis_Stag/Results/3_loop_o_gram/matrix_new/10K_norm/"
  #bedList = [bl1+"3B9_5-.txt:STAG1only_regions.txt", bl1+"3B9_5-.txt:STAG2only_regions.txt", bl1+"3B9_5plus.txt:STAG1only_regions.txt", bl1+"3B9_5plus.txt:STAG2only_regions.txt", bl1+"69_127-.txt:STAG1only_regions.txt", bl1+"69_127-.txt:STAG2only_regions.txt", bl1+"69_127plus.txt:STAG1only_regions.txt", bl1+"69_127plus.txt:STAG2only_regions.txt", bl2+"R_STAG1_RM_STAG1_intersection.bed", bl2+"R_STAG1_RM_STAG1_minusOnly.bed", bl2+"R_STAG1_RM_STAG1_plusOnly.bed"]

  bl2 = "/home/egg/Projects/Papantonis_Stag/Results/10_all_differential_loopograms/input/10K_norm/"
  bedList = [bl2+"R_STAG1_RM_STAG1_intersection.bed", bl2+"R_STAG1_RM_STAG1_minusOnly.bed", bl2+"R_STAG1_RM_STAG1_plusOnly.bed"]

  # Bed Loop
  for bedName in bedList:

    # Name
    ss = bamName.split("/")
    bname = ss[2]
    rn = bedName.split("/")[-1]
    if("R_STAG" in rn): rname = rn
    elif(rn == "3B9_5-.txt:STAG1only_regions.txt"): rname = "R_STAG1_RM_STAG2_minus"
    elif(rn == "3B9_5-.txt:STAG2only_regions.txt"): rname = "R_STAG2_RM_STAG2_minus"
    elif(rn == "3B9_5plus.txt:STAG1only_regions.txt"): rname = "R_STAG1_RM_STAG2_plus"
    elif(rn == "3B9_5plus.txt:STAG2only_regions.txt"): rname = "R_STAG2_RM_STAG2_plus"
    elif(rn == "69_127-.txt:STAG1only_regions.txt"): rname = "R_STAG1_RM_STAG1_minus"
    elif(rn == "69_127-.txt:STAG2only_regions.txt"): rname = "R_STAG2_RM_STAG1_minus"
    elif(rn == "69_127plus.txt:STAG1only_regions.txt"): rname = "R_STAG1_RM_STAG1_plus"
    elif(rn == "69_127plus.txt:STAG2only_regions.txt"): rname = "R_STAG2_RM_STAG1_plus"
    name = "__".join([bname, rname])
    if("DNase-seq" in bamName):
      extU = "5"
      extD = "5"
    else:
      extU = "0"
      extD = "200"

    # Input
    plotExt = "5000"
    signalExtU = extU
    signalExtD = extD
    bedFileName = bedName
    bamFileName = il+bamName+".bam"
    outputFileName = ol+name+".txt"

    # Execution
    print name
    command = "python 1_createTables.py "+" ".join([plotExt, signalExtU, signalExtD, bedFileName, bamFileName, outputFileName])
    os.system(command)


