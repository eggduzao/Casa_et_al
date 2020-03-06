
# Import
import os
import sys

# Region List
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/10_Losada_Meta_Plots_Genes/input/"
il = "/projects/ag-papan/eduardo/Wendt_Stag/Results/10_Losada_Meta_Plots_Genes/1_Meta_Gene_Tables/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/10_Losada_Meta_Plots_Genes/2_Meta_Gene_Plots/"
regionList = ["downregulated_by_STAG1_only", "downregulated_by_STAG1-STAG2", "downregulated_by_STAG2_only",
              "upregulated_by_STAG1_only", "upregulated_by_STAG1-STAG2", "upregulated_by_STAG2_only",
              "shared", "stag1_all", "stag1_only", "stag2_all", "stag2_only", "all_genes"]

# Open input file
inFileName = fl + "2_lgp.txt"
inFile = open(inFileName, "w")

# Region Loop
for region in regionList:

  # Signal List
  inLoc = il + region + "/"
  outLoc = ol + region + "/"
  command = "mkdir -p "+outLoc
  os.system(command)
  signalList = ["CTCF_MCF10A_CHIPSEQ", "SA1_MCF10A_CHIPSEQ", "SA2_MCF10A_CHIPSEQ", "SMC1_MCF10A_CHIPSEQ", "ZMYM2_MCF10A_CHIPSEQ"]

  # Signal Loop
  nameCounter = 1
  for signal in signalList:

    # Parameters
    signalName = signal.split("_")[0].upper()
    outFileName = "_".join([str(nameCounter), signalName])
    nameCounter += 1

    # Input
    label = signalName
    inputTableFileName = inLoc + signal + ".txt"
    outputAggrFileName = outLoc + outFileName + "_aggr.pdf"
    outputHeatFileName = outLoc + outFileName + "_heat.png"
    outputClusFileName = outLoc + outFileName + "_clus.txt"

    # Write
    inFile.write(" ".join([label, inputTableFileName, outputAggrFileName, outputHeatFileName, outputClusFileName])+"\n")

# Close file
inFile.close()


