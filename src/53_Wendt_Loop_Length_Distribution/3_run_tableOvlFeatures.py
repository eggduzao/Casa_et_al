
# Import
import os
import sys

# Folder List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/24_Wendt_Loop_Length_Distribution/input/"
ll = "/usr/users/egadegu/Projects/Wendt_Stag/Results/5_Loops/3_contact_files_filtered/"
rl1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/4_Peaks/"
rl2 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/rao/macs/"
rl3 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bed_files/"
tl = "/scratch/egadegu/CTC/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/24_Wendt_Loop_Length_Distribution/2_Categorized/"
loopNameList = ["STAG1-AUX", "STAG1+AUX", "STAG2-AUX", "STAG2+AUX"]
loopFileList = ["STAG1_minus_auxin_150Kb_5Mb" , "STAG1_plus_auxin_150Kb_5Mb" , "STAG2_minus_auxin_150Kb_5Mb" , "STAG2_plus_auxin_150Kb_5Mb"]
regionNameList = ["SMC3", "RAD21", "CTCF5", "CTCF7", "Shared", "Stag1", "Stag2"]
regionFileList = [rl1+"Smc3_peaks_filtered.narrowPeak", rl2+"CHIP001_RAD21_untreated_SE_summits_filtered.bed", rl2+"CHIP005_CTCF_untreated_PE_summits_filtered.bed", rl2+"CHIP007_CTCF_untreated_SE_summits_filtered.bed",
                  rl3+"shared.bed", rl3+"STAG1_only.bed", rl3+"STAG2_only.bed"]

# Opening file
inFileName = fl + "3_ctc.txt"
inFile = open(inFileName,"w")

# Input
intType = "1"
resolution = "25000"
halfExt = "100"
inputLoopNameList = ",".join(loopNameList)
inputLoopFileNameList = ",".join([ll + e + ".txt" for e in loopFileList])
inputPeaksNameList = ",".join(regionNameList)
inputPeaksFileNameList = ",".join(regionFileList)
temporaryLocation = tl
outputFileName = ol + "table.txt"

# Creating files
inFile.write(" ".join([intType, resolution, halfExt, inputLoopNameList, inputLoopFileNameList, inputPeaksNameList, inputPeaksFileNameList, temporaryLocation, outputFileName])+"\n")

# Closing file
inFile.close()


