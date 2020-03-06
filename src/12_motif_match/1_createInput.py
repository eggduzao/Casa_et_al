
# Import
import os
import sys
from glob import glob

# Motif List
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/12_motif_match/input/"
ml = "/projects/ag-papan/eduardo/Gusmao_DeNovo/Data/Gold_Standard/PFM_files/"
regionsFile = "/home/egusmao/rgtdata/hg19/binned_genome_100000.bed"
tl = "/scratch/eduardo/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/motif_match/e3/"
motifList = [ "H_TBP", "J_TBP", "H_CTCF", "J_CTCF", "J_TCF7", "J_TCF7L1", "J_TCF7L2", "H_TCF7L1", "H_TCF7L2", "H_TCF7", "h_ZNF274", "J_ZFX", "H_ZFX", "J_REST", "H_REST", "J_ZBTB33", "H_ZBTB33", "j_SP1", "h_SP1", "J_SRF", "H_SRF", "J_EGR1", "H_EGR1", "J_CEBPB", "H_CEBPB", "H_MAX", "J_MAX_MYC", "J_MAX", "H_FOSL1", "j_FOSL1_FOSL1_JUNB", "j_FOSL1_FOSL1_JUND", "J_FOSL1_JUN", "J_FOSL1", "h_TEAD4", "j_TEAD4", "J_ATF3", "H_ATF3", "H_USF1", "J_USF1", "J_YY1", "H_YY1", "J_ELF1", "H_ELF1"]

# Motif Loop
counter = 1
for motifName in motifList:

  # Parameters
  motifFileName = ml+motifName+".txt"
  regionsFileName = regionsFile
  tempLocation = tl+motifName+"/"
  outputFileName = ol+motifName+".bed"

  # Input File
  inFileName = il+str(counter)+".txt"
  inFile = open(inFileName, "w")
  inFile.write("\n".join([motifFileName, regionsFileName, tempLocation, outputFileName]))
  inFile.close()
  counter += 1


