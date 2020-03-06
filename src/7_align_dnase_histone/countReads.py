
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

# Input
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/rao/bam/"
bamFileNameList = ["CHIP001_RAD21_untreated_SE", "CHIP002_RAD21_treated_SE", "CHIP003_SMC1_untreated_PE", "CHIP004_SMC1_treated_PE", "CHIP005_CTCF_untreated_PE", "CHIP006_CTCF_treated_PE", "CHIP007_CTCF_untreated_SE", "CHIP008_CTCF_treated_SE", "CHIP009_H3K27Ac_untreated_SE", "CHIP010_H3K27Ac_treated_SE", "CHIP011_H3K4me3_untreated_SE", "CHIP012_H3K4me3_treated_SE", "CHIP013_H3K4me1_untreated_SE", "CHIP014_H3K4me1_treated_SE", "CHIP015_H3K36me3_untreated_SE", "CHIP016_H3K36me3_treated_SE", "CHIP017_H3K27me3_untreated_SE", "CHIP018_H3K27me3_treated_SE", "CHIP019_H3K9me3_untreated_SE", "CHIP020_H3K9me3_treated_SE", "CHIP021_H4K16Ac_untreated_SE", "CHIP022_H4K16Ac_treated_SE", "CHIP023_H3K79me2_untreated_SE", "CHIP024_H3K79me2_treated_SE", "CHIP025_H4K20me3_untreated_SE", "CHIP026_H4K20me3_treated_SE", "CHIP027_H2AZ_untreated_SE", "CHIP028_H2AZ_treated_SE", "CHIP029_NIPBL_untreated_SE", "CHIP030_NIPBL_treated_SE", "CHIP031_input_untreated_SE", "CHIP032_input_treated_SE", "CHIP033_input_untreated_PE", "CHIP034_input_treated_PE", "CHIP035_input_untreated_SE", "CHIP036_input_treated_SE", "CHIP037_input_untreated_SE", "CHIP038_input_treated_SE", "CHIP039_input_treated_SE"]
outputFileName = "count_list_rao.txt"

# Execution
outputFile = open(outputFileName,"w")
for bamFileName in bamFileNameList:
  tempFileName = "./a.txt"
  os.system("samtools view -c -F 260 "+il+bamFileName+".bam > "+tempFileName)
  tempFile = open(tempFileName,"r")
  value = tempFile.readline().strip()
  tempFile.close()
  outputFile.write("\t".join([bamFileName,value])+"\n")
  os.system("rm "+tempFileName)
outputFile.close()

"""
# Input
bl1 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/dnase/bam/"
bl2 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/histones/bam/"
bl3 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/tf/bam/"
bamFileNameList = [bl1+"HCT116_DNase-seq_UW.bam", bl2+"HCT116_ChIP-seq_H2AFZ_BROAD.bam", bl2+"HCT116_ChIP-seq_H3K4me1_USC.bam", bl2+"HCT116_ChIP-seq_H3K4me2_BROAD.bam", bl2+"HCT116_ChIP-seq_H3K4me3_UW.bam", bl2+"HCT116_ChIP-seq_H3K9ac_USC.bam", bl2+"HCT116_ChIP-seq_H3K9me2_BROAD.bam", bl2+"HCT116_ChIP-seq_H3K9me3_USC.bam", bl2+"HCT116_ChIP-seq_H3K20me1_BROAD.bam", bl2+"HCT116_ChIP-seq_H3K27ac_USC.bam", bl2+"HCT116_ChIP-seq_H3K27me3_BROAD.bam", bl2+"HCT116_ChIP-seq_H3K36me3_USC.bam", bl2+"HCT116_ChIP-seq_H3K79me2_BROAD.bam", bl3+"HCT116_ChIP-seq_ATF3_HAIB.bam", bl3+"HCT116_ChIP-seq_CBX3_HAIB.bam", bl3+"HCT116_ChIP-seq_CEBPB_HAIB.bam", bl3+"HCT116_ChIP-seq_CTCF_BROAD.bam", bl3+"HCT116_ChIP-seq_CTCF_HAIB.bam", bl3+"HCT116_ChIP-seq_CTCF_UW.bam", bl3+"HCT116_ChIP-seq_EGR1_HAIB.bam", bl3+"HCT116_ChIP-seq_ELF1_HAIB.bam", bl3+"HCT116_ChIP-seq_EZH2_BROAD.bam", bl3+"HCT116_ChIP-seq_EZH2phosphoT487_BROAD.bam", bl3+"HCT116_ChIP-seq_FOSL1_HAIB.bam", bl3+"HCT116_ChIP-seq_MAX_HAIB.bam", bl3+"HCT116_ChIP-seq_POL2RA_USC.bam", bl3+"HCT116_ChIP-seq_POLR2AphosphoS5_HAIB.bam", bl3+"HCT116_ChIP-seq_REST_HAIB.bam", bl3+"HCT116_ChIP-seq_SIN3A_HAIB.bam", bl3+"HCT116_ChIP-seq_SP1_HAIB.bam", bl3+"HCT116_ChIP-seq_SRF_HAIB.bam", bl3+"HCT116_ChIP-seq_TCF7L2_USC.bam", bl3+"HCT116_ChIP-seq_TEAD4_HAIB.bam", bl3+"HCT116_ChIP-seq_USF1_HAIB.bam", bl3+"HCT116_ChIP-seq_YY1_HAIB.bam", bl3+"HCT116_ChIP-seq_ZBTB33_HAIB.bam", bl3+"HCT116_ChIP-seq_ZNF274_USC.bam", bl3+"HCT116_ChIP-seq_JUND_HAIB.bam", bl3+"HCT116_ChIP-seq_RAD21_HAIB.bam", bl3+"HCT116_ChIP-seq_ZFX_SYDH.bam", bl3+"HCT116_ChIP-seq_ZFX_SYDH_PE_1_ENCFF811IJJ.bam", bl3+"HCT116_ChIP-seq_ZFX_SYDH_PE_2_ENCFF223VZW.bam"]
outputFileName = "count_list.txt"

# Execution
outputFile = open(outputFileName,"w")
for bamFileName in bamFileNameList:
  tempFileName = "./a.txt"
  os.system("samtools view -c -F 260 "+bamFileName+" > "+tempFileName)
  tempFile = open(tempFileName,"r")
  value = tempFile.readline().strip()
  tempFile.close()
  bamNameToWrite = bamFileName.split("/")[-1].split(".")[0]
  outputFile.write("\t".join([bamNameToWrite,value])+"\n")
  os.system("rm "+tempFileName)
outputFile.close()
"""

