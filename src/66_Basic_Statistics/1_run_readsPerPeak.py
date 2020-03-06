
# Import
import os
import sys

# Peak List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/37_Basic_Statistics/input/"
pl1 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bed_files/"
pl2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/4_Peaks/"
bl1 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bam_files/"
bl2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/2_Merged_Bam_Files/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/37_Basic_Statistics/1_ChIP/"
peakList = [pl1+"nonpredominant", pl1+"shared", pl1+"nonpredominant", pl1+"shared", pl1+"STAG1_full_peaks", pl1+"STAG1_only", pl1+"STAG1_predominant", pl1+"STAG2_full_peaks", pl1+"STAG2_only", pl1+"STAG2_predominant", pl2+"Smc3_peaks_filtered"]
bamList = [bl1+"STAG1_GFP.mrg-replicates.srt", bl1+"STAG1_GFP.mrg-replicates.srt", bl1+"STAG2_GFP.mrg-replicates.srt", bl1+"STAG2_GFP.mrg-replicates.srt", bl1+"STAG1_GFP.mrg-replicates.srt", bl1+"STAG1_GFP.mrg-replicates.srt", bl1+"STAG1_GFP.mrg-replicates.srt", bl1+"STAG2_GFP.mrg-replicates.srt", bl1+"STAG2_GFP.mrg-replicates.srt", bl1+"STAG2_GFP.mrg-replicates.srt", bl2+"Smc3"]

# Open File
inputFileName = fl + "1_rpp.txt"
inFile = open(inputFileName, "w")

# Peak Loop
for i in range(0, len(peakList)):

  # Parameters
  if("nonpredominant" in peakList[i] or "shared" in peakList[i]):
    name = peakList[i].split("/")[-1] + "_" + bamList[i].split("/")[-1].split("_")[0]
  else: name = peakList[i].split("/")[-1]
  if("Smc3_peaks_filtered" in peakList[i]): suff = ".narrowPeak"
  else: suff = ".bed"

  # Input
  regionFileName = peakList[i] + suff
  bamFileName = bamList[i] + ".bam"
  outputFileName = ol + name + ".txt"

  # Execution
  inFile.write(" ".join([regionFileName, bamFileName, outputFileName])+"\n")

inFile.close()


