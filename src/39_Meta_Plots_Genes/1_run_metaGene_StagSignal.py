
# Import
import os
import sys

# Region List
regionFile = "/home/egg/Projects/Wendt_Stag/Results/0_Definitive_Gene_Annotation/regions.bed"
aliasFile = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
rl = "/home/egg/Projects/Wendt_Stag/Data/expression/gene_lists/"
ol = "/home/egg/Projects/Wendt_Stag/Results/8_Meta_Plots_Genes/1_Meta_Gene_Tables/smoothed/"
regionList = ["all_downregulated_by_STAG1_and_specific", "all_downregulated_by_STAG2_and_specific",
              "all_upregulated_by_STAG2_and_specific", "downregulated_by_STAG2_only_and_specific",
              "upregulated_by_STAG2_only_and_specific", "downregulated_by_STAG1-STAG2_and_specific",
              "upregulated_by_STAG1-STAG2_and_specific", "downregulated_by_STAG1_only_and_specific",
              "all_upregulated_by_STAG1_and_specific", "upregulated_by_STAG1_only_and_specific"]

# Region Loop
for region in regionList:

  # Name
  outLoc = ol + "all_genes/"
  if(region != "."): outLoc = ol + region + "/"
  command = "mkdir -p "+outLoc
  os.system(command)

  # Signal List
  sl1 = "/home/egg/Projects/Wendt_Stag/Data/stag_bw_files/"
  sl2 = "/home/egg/Projects/Wendt_Stag/Data/rao/bam/"
  sl3 = "/home/egg/Projects/Wendt_Stag/Data/dnase/bam/"
  sl4 = "/home/egg/Projects/Wendt_Stag/Data/histones/bam/"
  sl5 = "/home/egg/Projects/Wendt_Stag/Data/tf/bam/"
  signalList = ["STAG1_GFP.mrg-replicates.srt.norm", "STAG2_GFP.mrg-replicates.srt.norm", "CHIP001_RAD21_untreated_SE", "CHIP003_SMC1_untreated_PE", "CHIP005_CTCF_untreated_PE", "CHIP007_CTCF_untreated_SE", "CHIP009_H3K27Ac_untreated_SE", "CHIP011_H3K4me3_untreated_SE", "CHIP013_H3K4me1_untreated_SE", "CHIP015_H3K36me3_untreated_SE", "CHIP017_H3K27me3_untreated_SE", "CHIP019_H3K9me3_untreated_SE", "CHIP021_H4K16Ac_untreated_SE", "CHIP023_H3K79me2_untreated_SE", "CHIP025_H4K20me3_untreated_SE", "CHIP027_H2AZ_untreated_SE", "CHIP029_NIPBL_untreated_SE", "HCT116_DNase-seq_UW", "HCT116_ChIP-seq_H2AFZ_BROAD", "HCT116_ChIP-seq_H3K4me1_USC", "HCT116_ChIP-seq_H3K4me2_BROAD", "HCT116_ChIP-seq_H3K4me3_UW", "HCT116_ChIP-seq_H3K9ac_USC", "HCT116_ChIP-seq_H3K9me2_BROAD", "HCT116_ChIP-seq_H3K9me3_USC", "HCT116_ChIP-seq_H3K27ac_USC", "HCT116_ChIP-seq_H3K27me3_BROAD", "HCT116_ChIP-seq_H3K36me3_USC", "HCT116_ChIP-seq_H3K79me2_BROAD", "HCT116_ChIP-seq_H4K20me1_BROAD", "HCT116_ChIP-seq_ATF3_HAIB", "HCT116_ChIP-seq_CBX3_HAIB", "HCT116_ChIP-seq_CEBPB_HAIB", "HCT116_ChIP-seq_CTCF_BROAD", "HCT116_ChIP-seq_CTCF_HAIB", "HCT116_ChIP-seq_CTCF_UW", "HCT116_ChIP-seq_EGR1_HAIB", "HCT116_ChIP-seq_ELF1_HAIB", "HCT116_ChIP-seq_EZH2_BROAD", "HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "HCT116_ChIP-seq_FOSL1_HAIB", "HCT116_ChIP-seq_JUND_HAIB", "HCT116_ChIP-seq_MAX_HAIB", "HCT116_ChIP-seq_POL2RA_USC", "HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "HCT116_ChIP-seq_RAD21_HAIB", "HCT116_ChIP-seq_REST_HAIB", "HCT116_ChIP-seq_SIN3A_HAIB", "HCT116_ChIP-seq_SP1_HAIB", "HCT116_ChIP-seq_SRF_HAIB", "HCT116_ChIP-seq_TCF7L2_USC", "HCT116_ChIP-seq_TEAD4_HAIB", "HCT116_ChIP-seq_USF1_HAIB", "HCT116_ChIP-seq_YY1_HAIB", "HCT116_ChIP-seq_ZBTB33_HAIB", "HCT116_ChIP-seq_ZFX_SYDH", "HCT116_ChIP-seq_ZNF274_USC"]
  countsList = ["20922666", "19917173", "24426210", "51504471", "69421213", "28035013", "12520401", "13930935", "25813227", "35677601", "26087129", "19096659", "21576671", "35721669", "23121517", "40007613", "30684077", "49075001", "27630623", "49435127", "114306062", "41097260", "53169949", "37660905", "55592854", "32818712", "25270406", "47271224", "24264417", "35483782", "48438364", "47462798", "42279850", "25415806", "19513799", "26856414", "31065624", "58707093", "36721190", "66380033", "36959453", "43802072", "50420217", "12346549", "35608075", "40122481", "40581337", "44235677", "46733429", "41374256", "13751559", "28998094", "77384084", "17350845", "25175554", "199586504", "48368885"]

  # Signal Loop
  nameCounter = 0
  for i in range(0,len(signalList)):

    # Parameters
    signal = signalList[i]
    ncounts = countsList[i]
    il = sl1
    ue = "0"
    de = "200"
    signalName = signal.split("_")[0].upper()
    signalType = "bw"
    if("CHIP0" in signal):
      il = sl2
      signalName = signal.split("_")[1].upper()
      signalType = "bam"
    elif("HCT116_DNase-seq" in signal):
      il = sl3
      ue = "5"
      de = "5"
      signalName = signal.split("_")[1].split("-")[0].upper()
      signalType = "bam"
    elif("HCT116_ChIP-seq_H" in signal):
      il = sl4
      signalName = signal.split("_")[2].upper()
      signalType = "bam"
    elif("HCT116" in signal):
      il = sl5
      signalName = signal.split("_")[2].upper()
      signalType = "bam"
    signalName = str(nameCounter) + "_" + signalName
    nameCounter += 1

    # Input
    downstream_extension = de
    upstream_extension = ue
    numberOfBins = "500"
    numberOfCounts = ncounts
    aliasFileName = aliasFile
    if(region == "."): geneListFileName = region
    else: geneListFileName = rl + region + ".txt"
    regionsBedFileName = regionFile
    signalFileType = signalType
    signalFileName = il + signal + "." + signalType
    outputGenesFileName = outLoc + signal + ".txt"

    # Write
    command = "python 1_metaGene_StagSignal.py "+" ".join([downstream_extension, upstream_extension, numberOfBins, numberOfCounts, aliasFileName, geneListFileName, regionsBedFileName, signalFileType, signalFileName, outputGenesFileName])+"\n"
    os.system(command)


