
# Import
import os
import sys

# ChIP/Motif Lists
fl = "./input_mwp/"
tl = "/scratch/eduardo/"
cl = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/tf/macs/"
ml = "/projects/ag-papan/eduardo/Gusmao_DeNovo/Data/Gold_Standard/Motif_Match/hg19/e3/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/tf/motifs/"
chipList = ["HCT116_ChIP-seq_ATF3_HAIB", "HCT116_ChIP-seq_CEBPB_HAIB", "HCT116_ChIP-seq_CTCF_HAIB", "HCT116_ChIP-seq_CTCF_UW", "HCT116_ChIP-seq_EGR1_HAIB", "HCT116_ChIP-seq_ELF1_HAIB", "HCT116_ChIP-seq_FOSL1_HAIB", "HCT116_ChIP-seq_MAX_HAIB", "HCT116_ChIP-seq_POL2RA_USC", "HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "HCT116_ChIP-seq_REST_HAIB", "HCT116_ChIP-seq_SP1_HAIB", "HCT116_ChIP-seq_TCF7L2_USC", "HCT116_ChIP-seq_TEAD4_HAIB", "HCT116_ChIP-seq_USF1_HAIB", "HCT116_ChIP-seq_YY1_HAIB", "HCT116_ChIP-seq_ZBTB33_HAIB", "HCT116_ChIP-seq_ZNF274_USC", "HCT116_ChIP-seq_SRF_HAIB", "HCT116_ChIP-seq_CTCF_BROAD"]
motifList = [["J_ATF3", "H_ATF3"], ["J_CEBPB", "H_CEBPB"], ["H_CTCF", "J_CTCF"], ["H_CTCF", "J_CTCF"], ["J_EGR1", "H_EGR1"], ["J_ELF1", "H_ELF1"], ["H_FOSL1", "j_FOSL1_FOSL1_JUNB", "j_FOSL1_FOSL1_JUND", "J_FOSL1_JUN", "J_FOSL1"], ["H_MAX", "J_MAX_MYC", "J_MAX"], ["H_TBP", "J_TBP"], ["H_TBP", "J_TBP"], ["J_REST", "H_REST"], ["j_SP1", "h_SP1"], ["J_TCF7L2", "H_TCF7L2", "H_TCF7", "J_TCF7", "J_TCF7L1", "H_TCF7L1"], ["h_TEAD4", "j_TEAD4"], ["H_USF1", "J_USF1"], ["J_YY1", "H_YY1"], ["J_ZBTB33", "H_ZBTB33"], ["h_ZNF274"], ["J_SRF", "H_SRF"], ["H_CTCF", "J_CTCF"]]

# ChIP/Motif Loop
counter = 1
for i in range(0,len(chipList)):

  # Parameters
  cName = chipList[i]
  mName = ",".join([ml+e+".bed" for e in motifList[i]])

  # Input
  motifFileNameList = mName
  chromSizesFileName = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
  summitFileName = cl+cName+"_summits_filter.bed"
  tempLocation = tl+cName+"_MWP/"
  outputFileName = ol+cName+".bed"

  # Execution
  inputFileName = fl+str(counter)+".txt"
  inputFile = open(inputFileName, "w")
  inputFile.write("\n".join([motifFileNameList, chromSizesFileName, summitFileName, tempLocation, outputFileName]))
  inputFile.close()
  counter += 1


