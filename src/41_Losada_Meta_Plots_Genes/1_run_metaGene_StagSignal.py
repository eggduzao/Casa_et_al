
# Import
import os
import sys

# Region List
regionFile = "/projects/ag-papan/eduardo/Wendt_Stag/Results/0_Definitive_Gene_Annotation/regions.bed"
aliasFile = "/home/egusmao/rgtdata/hg19/alias_human_booster.txt"
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/10_Losada_Meta_Plots_Genes/input/"
rl = "/projects/ag-papan/eduardo/Wendt_Stag/Data/losada/expression/gene_lists/"
sl = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/2_Merged_Bam_Files/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/10_Losada_Meta_Plots_Genes/1_Meta_Gene_Tables/"
regionList = ["downregulated_by_STAG1_only", "downregulated_by_STAG1-STAG2", "downregulated_by_STAG2_only",
              "upregulated_by_STAG1_only", "upregulated_by_STAG1-STAG2", "upregulated_by_STAG2_only",
              "shared", "stag1_all", "stag1_only", "stag2_all", "stag2_only", "."]

# Write
inFileName = fl + "1_lgs.txt"
inFile = open(inFileName, "w")

# Region Loop
for region in regionList:

  # Name
  outLoc = ol + "all_genes/"
  if(region != "."): outLoc = ol + region + "/"
  command = "mkdir -p "+outLoc
  os.system(command)

  # Signal List
  signalList = ["CTCF_MCF10A_CHIPSEQ", "SA1_MCF10A_CHIPSEQ", "SA2_MCF10A_CHIPSEQ", "SMC1_MCF10A_CHIPSEQ", "ZMYM2_MCF10A_CHIPSEQ"]
  countsList = ["34917088", "35600693", "38385254", "45107605", "17630609"]

  # Signal Loop
  for i in range(0,len(signalList)):
  
    # Parameters
    signal = signalList[i]
    ncounts = countsList[i]

    # Input
    downstream_extension = "200"
    upstream_extension = "0"
    numberOfBins = "100"
    numberOfCounts = ncounts
    aliasFileName = aliasFile
    if(region == "."): geneListFileName = region
    else: geneListFileName = rl + region + ".txt"
    regionsBedFileName = regionFile
    signalFileType = "bam"
    signalFileName = sl + signal + ".bam"
    outputGenesFileName = outLoc + signal + ".txt"

    # Write
    inFile.write(" ".join([downstream_extension, upstream_extension, numberOfBins, numberOfCounts, aliasFileName, geneListFileName, regionsBedFileName, signalFileType, signalFileName, outputGenesFileName])+"\n")

# Closing file
inFile.close()


