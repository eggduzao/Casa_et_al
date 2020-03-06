
# Import
import os
import sys

# Region List
il = "/home/egg/Projects/Wendt_Stag/Results/8_Meta_Plots_Genes/1_Meta_Gene_Tables/smoothed/"
ol = "/home/egg/Projects/Wendt_Stag/Results/8_Meta_Plots_Genes/2_Meta_Gene_Plots/smoothed/"
regionList = ["all_downregulated_by_STAG1_and_specific", "all_downregulated_by_STAG2_and_specific",
              "all_upregulated_by_STAG2_and_specific", "downregulated_by_STAG2_only_and_specific",
              "upregulated_by_STAG2_only_and_specific", "downregulated_by_STAG1-STAG2_and_specific",
              "upregulated_by_STAG1-STAG2_and_specific", "downregulated_by_STAG1_only_and_specific",
              "all_upregulated_by_STAG1_and_specific", "upregulated_by_STAG1_only_and_specific"]

# Region Loop
for region in regionList:

  # Signal List
  inLoc = il + region + "/"
  outLoc = ol + region + "/"
  command = "mkdir -p "+outLoc
  os.system(command)
  signalList = ["STAG1_GFP.mrg-replicates.srt.norm", "STAG2_GFP.mrg-replicates.srt.norm", "CHIP001_RAD21_untreated_SE", "CHIP003_SMC1_untreated_PE", "CHIP007_CTCF_untreated_SE", "CHIP009_H3K27Ac_untreated_SE", "CHIP011_H3K4me3_untreated_SE", "CHIP013_H3K4me1_untreated_SE", "CHIP015_H3K36me3_untreated_SE", "CHIP017_H3K27me3_untreated_SE", "CHIP019_H3K9me3_untreated_SE", "CHIP021_H4K16Ac_untreated_SE", "CHIP023_H3K79me2_untreated_SE", "CHIP025_H4K20me3_untreated_SE", "CHIP027_H2AZ_untreated_SE", "CHIP029_NIPBL_untreated_SE", "HCT116_DNase-seq_UW", "HCT116_ChIP-seq_H2AFZ_BROAD", "HCT116_ChIP-seq_H3K4me1_USC", "HCT116_ChIP-seq_H3K4me2_BROAD", "HCT116_ChIP-seq_H3K4me3_UW", "HCT116_ChIP-seq_H3K9ac_USC", "HCT116_ChIP-seq_H3K9me2_BROAD", "HCT116_ChIP-seq_H3K9me3_USC", "HCT116_ChIP-seq_H3K27ac_USC", "HCT116_ChIP-seq_H3K27me3_BROAD", "HCT116_ChIP-seq_H3K36me3_USC", "HCT116_ChIP-seq_H3K79me2_BROAD", "HCT116_ChIP-seq_H4K20me1_BROAD", "HCT116_ChIP-seq_ATF3_HAIB", "HCT116_ChIP-seq_CBX3_HAIB", "HCT116_ChIP-seq_CEBPB_HAIB", "HCT116_ChIP-seq_CTCF_BROAD", "HCT116_ChIP-seq_CTCF_HAIB", "HCT116_ChIP-seq_CTCF_UW", "HCT116_ChIP-seq_EGR1_HAIB", "HCT116_ChIP-seq_ELF1_HAIB", "HCT116_ChIP-seq_EZH2_BROAD", "HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "HCT116_ChIP-seq_FOSL1_HAIB", "HCT116_ChIP-seq_JUND_HAIB", "HCT116_ChIP-seq_MAX_HAIB", "HCT116_ChIP-seq_POL2RA_USC", "HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "HCT116_ChIP-seq_RAD21_HAIB", "HCT116_ChIP-seq_REST_HAIB", "HCT116_ChIP-seq_SIN3A_HAIB", "HCT116_ChIP-seq_SP1_HAIB", "HCT116_ChIP-seq_SRF_HAIB", "HCT116_ChIP-seq_TCF7L2_USC", "HCT116_ChIP-seq_TEAD4_HAIB", "HCT116_ChIP-seq_USF1_HAIB", "HCT116_ChIP-seq_YY1_HAIB", "HCT116_ChIP-seq_ZBTB33_HAIB", "HCT116_ChIP-seq_ZNF274_USC"]
# "HCT116_ChIP-seq_ZFX_SYDH",

  # Signal Loop
  nameCounter = 1
  for signal in signalList:

    # Parameters
    if("GFP.mrg-replicates" in signal):
      signalName = signal.split("_")[0].upper()
      outFileName = "_".join([str(nameCounter), signalName])
    elif("CHIP0" in signal):
      signalName = "_".join([signal.split("_")[1].upper(), "Rao"])
      outFileName = "_".join([str(nameCounter), signal.split("_")[1].upper(), "RAO"])
    elif("HCT116_DNase-seq" in signal):
      signalName = "DNase_UW"
      outFileName = "_".join([str(nameCounter), signalName, "UW"])
    elif("HCT116" in signal):
      signalName = "_".join([signal.split("_")[2].upper(), signal.split("_")[-1].upper()])
      outFileName = "_".join([str(nameCounter), signal.split("_")[2].upper(), signal.split("_")[-1].upper()])
    nameCounter += 1

    # Input
    label = signalName
    inputTableFileName = inLoc + signal + ".txt"
    outputAggrFileName = outLoc + outFileName + "_aggr.pdf"
    outputHeatFileName = outLoc + outFileName + "_heat.png"
    outputClusFileName = outLoc + outFileName + "_clus.txt"

    # Write
    command = "R CMD BATCH '--args '"+label+"' '"+inputTableFileName+"' '"+outputAggrFileName+"' '"+outputHeatFileName+"' '"+outputClusFileName+" 2_metaGenePlot.R 2_metaGenePlot.Rout"
    os.system(command)


