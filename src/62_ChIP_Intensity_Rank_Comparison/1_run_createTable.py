
# Import
import os
import sys

# Stag List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/33_ChIP_Intensity_Rank_Comparison/input/"
wpl = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bed_files/"
wbl = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bam_files/"
lpl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/4_Peaks/"
lbl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/2_Merged_Bam_Files/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/33_ChIP_Intensity_Rank_Comparison/1_Tables/"
peakList = ["nonpredominant", "nonpredominant", "shared", "shared", "STAG1_full_peaks", "STAG1_only", "STAG1_predominant", "STAG2_full_peaks", "STAG2_only", "STAG2_predominant",
            "SA1_MCF10A_ChIP-seq", "SA1_MCF10A_ChIP-seq_siSA1", "SA1_MCF10A_ChIP-seq_siSA2", "SA2_MCF10A_ChIP-seq", "SA2_MCF10A_ChIP-seq_siSA1", "SA2_MCF10A_ChIP-seq_siSA2"]
bamList = ["STAG1", "STAG2", "STAG1", "STAG2", "STAG1", "STAG1", "STAG1", "STAG2", "STAG2", "STAG2",
             "SA1_MCF10A_ChIP-seq", "SA1_MCF10A_ChIP-seq_siSA1", "SA1_MCF10A_ChIP-seq_siSA2", "SA2_MCF10A_ChIP-seq", "SA2_MCF10A_ChIP-seq_siSA1", "SA2_MCF10A_ChIP-seq_siSA2"]
totalCountList = ["27008128", "25856287", "27008128", "25856287", "27008128", "27008128", "27008128", "25856287", "25856287", "25856287",
              "35600679", "40284495", "40739155", "38385265", "40884461", "45607748"]

# Open File
inputFileName = fl + "1_cic.txt"
inFile = open(inputFileName, "w")

# Stag Loop
for i in range(0,len(peakList)):

  # Paramters
  if("STAG" in bamList[i]):
    pl = wpl
    bl = wbl
    suffp = ".bed"
    suffb = "_GFP.mrg-replicates.srt.bam"
  else:
    pl = lpl
    bl = lbl
    suffp = "_peaks_filtered.narrowPeak"
    suffb = ".bam"

  # Input
  halfExt = "500"
  totalCount = totalCountList[i]
  regionFileName = pl + peakList[i] + suffp
  bamFileName = bl + bamList[i] + suffb
  outputFileName = ol + peakList[i] + ".txt"

  # Execution
  inFile.write(" ".join([halfExt, totalCount, regionFileName, bamFileName, outputFileName])+"\n")

inFile.close()


