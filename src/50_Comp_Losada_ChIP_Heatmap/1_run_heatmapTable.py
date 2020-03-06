
# Import
import os
import sys

# Folder List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/21_Comp_Losada_ChIP_Heatmap/input/"
bedL1 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bed_files/"
bedL2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/4_Peaks/"
bamL1 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bam_files/"
bamL2 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/2_Merged_Bam_Files/"
tl = "/scratch/egadegu/WLH/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/21_Comp_Losada_ChIP_Heatmap/1_Table/"
bedPairList = [["STAG1_full_peaks", "SA1_MCF10A_ChIP-seq_peaks_filtered"], ["STAG1_full_peaks", "SA2_MCF10A_ChIP-seq_peaks_filtered"], ["STAG2_full_peaks", "SA1_MCF10A_ChIP-seq_peaks_filtered"], ["STAG2_full_peaks", "SA2_MCF10A_ChIP-seq_peaks_filtered"]]
bamPairList = [["STAG1_GFP.mrg-replicates.srt", "SA1_MCF10A_ChIP-seq"], ["STAG1_GFP.mrg-replicates.srt", "SA2_MCF10A_ChIP-seq"], ["STAG2_GFP.mrg-replicates.srt", "SA1_MCF10A_ChIP-seq"], ["STAG2_GFP.mrg-replicates.srt", "SA2_MCF10A_ChIP-seq"]]
readPairList = [["27008128", "35600679"], ["27008128", "38385265"], ["25856287", "35600679"], ["25856287", "38385265"]]

# Opening file
inFileName = fl + "1_wlh.txt"
inFile = open(inFileName,"w")

# Fasta Loop
for i in range(0, len(bedPairList)):

  # Parameters
  bedPair = bedPairList[i]
  bamPair = bamPairList[i]
  readPair = readPairList[i]
  if(bedPair[0].split("_")[0] == "STAG1"): saW = "SA1"
  else: saW = "SA2"
  saL = bedPair[1].split("_")[0]
  name = "W_" + saW + "_L_" + saL

  # Input
  numberOfBins = "250"
  totalExtension = "2500"
  totalReadsW = readPair[0]
  totalReadsL = readPair[1]
  bedFileNameW = bedL1 + bedPair[0] + ".bed"
  bedFileNameL = bedL2 + bedPair[1] + ".narrowPeak"
  bamFileNameW = bamL1 + bamPair[0] + ".bam"
  bamFileNameL = bamL2 + bamPair[1] + ".bam"
  tempLocation = tl + name + "/"
  outputFileName = ol + name + "_2.txt"

  # Creating files
  inFile.write(" ".join([numberOfBins, totalExtension, totalReadsW, totalReadsL, bedFileNameW, bedFileNameL, bamFileNameW, bamFileNameL, tempLocation, outputFileName])+"\n")

# Closing file
inFile.close()

"""
# Input
numberOfBins = "50"
totalExtension = "2500"
bedFileNameW = bedL1 + bedPairList[0][0] + ".bed"
bedFileNameL = bedL2 + bedPairList[0][1] + ".narrowPeak"
bamFileNameW = bamL1 + bamPairList[0][0] + ".bam"
bamFileNameL = bamL2 + bamPairList[0][1] + ".bam"
tempLocation = "./test/"
outputFileName = "./res/a.txt"

os.system("python /usr/users/egadegu/Projects/Wendt_Stag/Code/21_Comp_Losada_ChIP_Heatmap/1_heatmapTable.py " + " ".join([numberOfBins, totalExtension, bedFileNameW, bedFileNameL, bamFileNameW, bamFileNameL, tempLocation, outputFileName]))
"""


