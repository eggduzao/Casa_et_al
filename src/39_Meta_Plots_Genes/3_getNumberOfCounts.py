
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Signal List
sl1 = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/stag_bam_files/"
sl2 = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/rao/bam/"
sl3 = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/dnase/bam/"
sl4 = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/histones/bam/"
sl5 = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/tf/bam/"
bamFileNameList = ["STAG1_GFP.mrg-replicates.srt", "STAG2_GFP.mrg-replicates.srt", "CHIP001_RAD21_untreated_SE", "CHIP003_SMC1_untreated_PE", "CHIP005_CTCF_untreated_PE", "CHIP007_CTCF_untreated_SE", "CHIP009_H3K27Ac_untreated_SE", "CHIP011_H3K4me3_untreated_SE", "CHIP013_H3K4me1_untreated_SE", "CHIP015_H3K36me3_untreated_SE", "CHIP017_H3K27me3_untreated_SE", "CHIP019_H3K9me3_untreated_SE", "CHIP021_H4K16Ac_untreated_SE", "CHIP023_H3K79me2_untreated_SE", "CHIP025_H4K20me3_untreated_SE", "CHIP027_H2AZ_untreated_SE", "CHIP029_NIPBL_untreated_SE", "HCT116_DNase-seq_UW", "HCT116_ChIP-seq_H2AFZ_BROAD", "HCT116_ChIP-seq_H3K4me1_USC", "HCT116_ChIP-seq_H3K4me2_BROAD", "HCT116_ChIP-seq_H3K4me3_UW", "HCT116_ChIP-seq_H3K9ac_USC", "HCT116_ChIP-seq_H3K9me2_BROAD", "HCT116_ChIP-seq_H3K9me3_USC", "HCT116_ChIP-seq_H3K27ac_USC", "HCT116_ChIP-seq_H3K27me3_BROAD", "HCT116_ChIP-seq_H3K36me3_USC", "HCT116_ChIP-seq_H3K79me2_BROAD", "HCT116_ChIP-seq_H4K20me1_BROAD", "HCT116_ChIP-seq_ATF3_HAIB", "HCT116_ChIP-seq_CBX3_HAIB", "HCT116_ChIP-seq_CEBPB_HAIB", "HCT116_ChIP-seq_CTCF_BROAD", "HCT116_ChIP-seq_CTCF_HAIB", "HCT116_ChIP-seq_CTCF_UW", "HCT116_ChIP-seq_EGR1_HAIB", "HCT116_ChIP-seq_ELF1_HAIB", "HCT116_ChIP-seq_EZH2_BROAD", "HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "HCT116_ChIP-seq_FOSL1_HAIB", "HCT116_ChIP-seq_JUND_HAIB", "HCT116_ChIP-seq_MAX_HAIB", "HCT116_ChIP-seq_POL2RA_USC", "HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "HCT116_ChIP-seq_RAD21_HAIB", "HCT116_ChIP-seq_REST_HAIB", "HCT116_ChIP-seq_SIN3A_HAIB", "HCT116_ChIP-seq_SP1_HAIB", "HCT116_ChIP-seq_SRF_HAIB", "HCT116_ChIP-seq_TCF7L2_USC", "HCT116_ChIP-seq_TEAD4_HAIB", "HCT116_ChIP-seq_USF1_HAIB", "HCT116_ChIP-seq_YY1_HAIB", "HCT116_ChIP-seq_ZBTB33_HAIB", "HCT116_ChIP-seq_ZFX_SYDH", "HCT116_ChIP-seq_ZNF274_USC"]
outputFileName = "./number_of_valid_reads.txt"

# Execution
outputFile = open(outputFileName,"w")
for bamName in bamFileNameList:
  il = sl1
  if("CHIP0" in bamName): il = sl2
  elif("HCT116_DNase-seq" in bamName): il = sl3
  elif("HCT116_ChIP-seq_H" in bamName): il = sl4
  elif("HCT116" in bamName): il = sl5
  tempFileName = "./a.txt"
  os.system("samtools view -c -F 260 "+ il + bamName + ".bam > " + tempFileName)
  tempFile = open(tempFileName,"rU")
  value = tempFile.readline().strip()
  tempFile.close()
  outputFile.write("\t".join([bamName,value])+"\n")
  os.system("rm "+tempFileName)
outputFile.close()


