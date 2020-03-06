
# Import
import os
import sys

# Stag List
#counter = 1
#fl = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/15_big_statistics_table/input_stab/"
sl1 = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/input/10K_norm/"
sl2 = "/home/egg/Projects/Wendt_Stag/Data/stag_bed_files/"
bl = "/home/egg/Projects/Wendt_Stag/Data/rao/bam/"
bl1 = "/home/egg/Projects/Wendt_Stag/Data/dnase/bam/"
bl2 = "/home/egg/Projects/Wendt_Stag/Data/histones/bam/"
bl3 = "/home/egg/Projects/Wendt_Stag/Data/tf/bam/"
bl4 = "/home/egg/Projects/Wendt_Stag/Data/rao/bam/"
cl = "/home/egg/Projects/Wendt_Stag/Data/control/bam/"
cl1 = "/home/egg/Projects/Wendt_Stag/Data/dnase/macs/"
cl2 = "/home/egg/Projects/Wendt_Stag/Data/histones/macs/"
cl3 = "/home/egg/Projects/Wendt_Stag/Data/tf/macs/"
cl4 = "/home/egg/Projects/Wendt_Stag/Data/rao/macs/"
el = "/home/egg/Projects/Wendt_Stag/Data/expression/"
ml = "/home/egg/Projects/Wendt_Stag/Data/tf/motifs/"
ol = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/"
#stagFileList = ["R_STAG1_RM_STAG1_intersection", "R_STAG1_RM_STAG1_minusOnly", "R_STAG1_RM_STAG1_plusOnly", "R_STAG2_RM_STAG2_intersection", "R_STAG2_RM_STAG2_minusOnly", "R_STAG2_RM_STAG2_plusOnly", "STAG1_only", "STAG2_only"]
stagFileList = ["R_STAG2_RM_STAG2_intersection", "R_STAG2_RM_STAG2_minusOnly", "R_STAG2_RM_STAG2_plusOnly"]

# Stag Loop
for stagFile in stagFileList:

  if("R_" in stagFile):
    il = sl1
    fileIL = "Y"
  else:
    il = sl2
    fileIL = "N"

  # Input
  fileIsLoop = fileIL
  stagFileName = il + stagFile + ".bed"
  aliasFileName = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
  chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
  genomeFileName = "/home/egg/Data/Genomes/hg19/hg19.fa"
  regionsFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/0_genomicRegions/hg19_regions.bam"
  ensemblDictFileName = "/home/egg/Projects/Wendt_Stag/Data/gene_annotation/hg19_ensembl_dictionary_filtered.txt"
  enhancersFileName = "/home/egg/Projects/Wendt_Stag/Data/gene_annotation/HCT116_enhancers_filter.bam"
  chrommHmmFileName = "NOTAVAILABLE"
  expressionLabelList = "STAG1,STAG2"
  expressionFileNameList = el+"STAG1.csv,"+el+"STAG2.csv"

  signalCountList1 = ",".join(["49075001", "27630623", "49435127", "114306062", "41097260", "53169949", "37660905", "55592854", "35483782", "32818712", "25270406", "47271224", "24264417", "48438364", "47462798", "42279850", "25415806", "19513799", "26856414", "31065624", "58707093", "36721190", "66380033", "36959453", "50420217", "12346549", "35608075", "40581337", "44235677", "46733429", "41374256", "13751559", "28998094", "77384084", "17350845", "25175554", "48368885", "40122481"])
  signalCountList2 = ",".join(["24426210", "51504471", "69421213", "28035013"])
  signalCountList = ",".join([signalCountList1, signalCountList2])

  signalFileNameList1 = bl1+"HCT116_DNase-seq_UW.bam,"+bl2+"HCT116_ChIP-seq_H2AFZ_BROAD.bam,"+bl2+"HCT116_ChIP-seq_H3K4me1_USC.bam,"+bl2+"HCT116_ChIP-seq_H3K4me2_BROAD.bam,"+bl2+"HCT116_ChIP-seq_H3K4me3_UW.bam,"+bl2+"HCT116_ChIP-seq_H3K9ac_USC.bam,"+bl2+"HCT116_ChIP-seq_H3K9me2_BROAD.bam,"+bl2+"HCT116_ChIP-seq_H3K9me3_USC.bam,"+bl2+"HCT116_ChIP-seq_H4K20me1_BROAD.bam,"+bl2+"HCT116_ChIP-seq_H3K27ac_USC.bam,"+bl2+"HCT116_ChIP-seq_H3K27me3_BROAD.bam,"+bl2+"HCT116_ChIP-seq_H3K36me3_USC.bam,"+bl2+"HCT116_ChIP-seq_H3K79me2_BROAD.bam,"+bl3+"HCT116_ChIP-seq_ATF3_HAIB.bam,"+bl3+"HCT116_ChIP-seq_CBX3_HAIB.bam,"+bl3+"HCT116_ChIP-seq_CEBPB_HAIB.bam,"+bl3+"HCT116_ChIP-seq_CTCF_BROAD.bam,"+bl3+"HCT116_ChIP-seq_CTCF_HAIB.bam,"+bl3+"HCT116_ChIP-seq_CTCF_UW.bam,"+bl3+"HCT116_ChIP-seq_EGR1_HAIB.bam,"+bl3+"HCT116_ChIP-seq_ELF1_HAIB.bam,"+bl3+"HCT116_ChIP-seq_EZH2_BROAD.bam,"+bl3+"HCT116_ChIP-seq_EZH2phosphoT487_BROAD.bam,"+bl3+"HCT116_ChIP-seq_FOSL1_HAIB.bam,"+bl3+"HCT116_ChIP-seq_MAX_HAIB.bam,"+bl3+"HCT116_ChIP-seq_POL2RA_USC.bam,"+bl3+"HCT116_ChIP-seq_POLR2AphosphoS5_HAIB.bam,"+bl3+"HCT116_ChIP-seq_REST_HAIB.bam,"+bl3+"HCT116_ChIP-seq_SIN3A_HAIB.bam,"+bl3+"HCT116_ChIP-seq_SP1_HAIB.bam,"+bl3+"HCT116_ChIP-seq_SRF_HAIB.bam,"+bl3+"HCT116_ChIP-seq_TCF7L2_USC.bam,"+bl3+"HCT116_ChIP-seq_TEAD4_HAIB.bam,"+bl3+"HCT116_ChIP-seq_USF1_HAIB.bam,"+bl3+"HCT116_ChIP-seq_YY1_HAIB.bam,"+bl3+"HCT116_ChIP-seq_ZBTB33_HAIB.bam,"+bl3+"HCT116_ChIP-seq_ZNF274_USC.bam,"+bl3+"HCT116_ChIP-seq_RAD21_HAIB.bam"
  signalFileNameList2 = bl4+"CHIP001_RAD21_untreated_SE.bam,"+bl4+"CHIP003_SMC1_untreated_PE.bam,"+bl4+"CHIP005_CTCF_untreated_PE.bam,"+bl4+"CHIP007_CTCF_untreated_SE.bam"
  signalFileNameList = ",".join([signalFileNameList1, signalFileNameList2])

  controlCountList = ",".join(["8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745", "8621745"])

  controlFileNameList = cl+"HCT116_ChIP-seq_Control_UW.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_USC.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_UW.bam,"+cl+"HCT116_ChIP-seq_Control_USC.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_USC.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_USC.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_USC.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_UW.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_BROAD.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_USC.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_USC.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+cl+"HCT116_ChIP-seq_Control_USC.bam,"+cl+"HCT116_ChIP-seq_Control_HAIB.bam,"+bl+"CHIP031_input_untreated_SE.bam,"+bl+"CHIP033_input_untreated_PE.bam,"+bl+"CHIP033_input_untreated_PE.bam,"+bl+"CHIP031_input_untreated_SE.bam"

  peakFileNameList = cl1+"HCT116_DNase-seq_UW_peaks_filter500.narrowPeak,"+cl2+"HCT116_ChIP-seq_H2AFZ_BROAD_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K4me1_USC_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K4me2_BROAD_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K4me3_UW_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K9ac_USC_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K9me2_BROAD_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K9me3_USC_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H4K20me1_BROAD_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K27ac_USC_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K27me3_BROAD_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K36me3_USC_peaks_filter.narrowPeak,"+cl2+"HCT116_ChIP-seq_H3K79me2_BROAD_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_ATF3_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_CBX3_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_CEBPB_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_CTCF_BROAD_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_CTCF_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_CTCF_UW_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_EGR1_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_ELF1_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_EZH2_BROAD_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_EZH2phosphoT487_BROAD_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_FOSL1_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_MAX_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_POL2RA_USC_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_POLR2AphosphoS5_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_REST_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_SIN3A_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_SP1_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_SRF_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_TCF7L2_USC_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_TEAD4_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_USF1_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_YY1_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_ZBTB33_HAIB_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_ZNF274_USC_peaks_filter.narrowPeak,"+cl3+"HCT116_ChIP-seq_RAD21_HAIB_peaks_filtered.narrowPeak,"+cl4+"CHIP001_RAD21_untreated_SE_peaks_filtered.narrowPeak,"+cl4+"CHIP003_SMC1_untreated_PE_peaks_filtered.narrowPeak,"+cl4+"CHIP005_CTCF_untreated_PE_peaks_filtered.narrowPeak,"+cl4+"CHIP007_CTCF_untreated_SE_peaks_filtered.narrowPeak"

  signalLabelList = ",".join(["_".join(e.split("/")[-1].split(".")[0].split("_")[-2:]) for e in signalFileNameList1.split(",")]) + "," + ",".join(["_".join(e.split("/")[-1].split(".")[0].split("_")[:2]) for e in signalFileNameList2.split(",")])
  motifFileNameList = ml+"HCT116_ChIP-seq_ATF3_HAIB.bam,"+ml+"HCT116_ChIP-seq_CEBPB_HAIB.bam,"+ml+"HCT116_ChIP-seq_CTCF_BROAD.bam,"+ml+"HCT116_ChIP-seq_CTCF_HAIB.bam,"+ml+"HCT116_ChIP-seq_CTCF_UW.bam,"+ml+"HCT116_ChIP-seq_EGR1_HAIB.bam,"+ml+"HCT116_ChIP-seq_ELF1_HAIB.bam,"+ml+"HCT116_ChIP-seq_FOSL1_HAIB.bam,"+ml+"HCT116_ChIP-seq_MAX_HAIB.bam,"+ml+"HCT116_ChIP-seq_POL2RA_USC.bam,"+ml+"HCT116_ChIP-seq_POLR2AphosphoS5_HAIB.bam,"+ml+"HCT116_ChIP-seq_REST_HAIB.bam,"+ml+"HCT116_ChIP-seq_SP1_HAIB.bam,"+ml+"HCT116_ChIP-seq_SRF_HAIB.bam,"+ml+"HCT116_ChIP-seq_TCF7L2_USC.bam,"+ml+"HCT116_ChIP-seq_TEAD4_HAIB.bam,"+ml+"HCT116_ChIP-seq_USF1_HAIB.bam,"+ml+"HCT116_ChIP-seq_YY1_HAIB.bam,"+ml+"HCT116_ChIP-seq_ZBTB33_HAIB.bam,"+ml+"HCT116_ChIP-seq_ZNF274_USC.bam,"+ml+"HCT116_ChIP-seq_ZNF274_USC.bam,"+ml+"HCT116_ChIP-seq_ZNF274_USC.bam,"+ml+"HCT116_ChIP-seq_ZNF274_USC.bam,"+ml+"HCT116_ChIP-seq_ZNF274_USC.bam,"+ml+"HCT116_ChIP-seq_ZNF274_USC.bam"
  motifLabelList = ",".join(["_".join(e.split("/")[-1].split(".")[0].split("_")[-2:]) for e in motifFileNameList.split(",")])
  tempLocation = "./TEMP/"
  outputFileName1 = ol+stagFile+"_genome_table.csv"
  outputFileName2 = ol+stagFile+"_signal_table.csv"
  outputFileName3 = ol+stagFile+"_motif_table.csv"
  outputFileName4 = ol+stagFile+"_peak_signal_table.csv"

  # Creating input file

  #inputFileName = fl+str(counter)+".txt"
  #inputFile = open(inputFileName, "w")
  #inputFile.write("\n".join([fileIsLoop, stagFileName, aliasFileName, chromSizesFileName, genomeFileName, regionsFileName, ensemblDictFileName, enhancersFileName, chrommHmmFileName, expressionLabelList, expressionFileNameList, signalLabelList, signalCountList, si#gnalFileNameList, controlCountList, controlFileNameList, peakFileNameList, motifLabelList, motifFileNameList, tempLocation, outputFileName1, outputFileName2, outputFileName3, outputFileName4]))
  #inputFile.close()
  #counter += 1

  command = "python 1_createTable.py "+" ".join([fileIsLoop, stagFileName, aliasFileName, chromSizesFileName, genomeFileName, regionsFileName, ensemblDictFileName, enhancersFileName, chrommHmmFileName, expressionLabelList, expressionFileNameList, signalLabelList, signalCountList, signalFileNameList, controlCountList, controlFileNameList, peakFileNameList, motifLabelList, motifFileNameList, tempLocation, outputFileName1, outputFileName2, outputFileName3, outputFileName4])
  os.system(command)


