
# Import
import os
import sys

# Stag List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/32_AB_TAD_Enrichments/input/"
gl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/0_Definitive_Gene_Annotation/"
sl = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bed_files/"
cl = "/usr/users/egadegu/Projects/Wendt_Stag/Data/rao/macs/"
tl = "/usr/users/egadegu/Projects/Wendt_Stag/Results/4_TADs/1_TADs_GMAP/"
al = "/usr/users/egadegu/Projects/Wendt_Stag/Results/3_AB_Compartments/2_compartments/"
#ml = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_matrix_files/25K_norm/" # Normal
ml = "/usr/users/egadegu/Projects/Wendt_Stag/Results/30_Subtraction_Maps/1_Bin_Std_Matrices/25K_norm/" # Binned
tempLoc = "/scratch/egadegu/ATEB/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/32_AB_TAD_Enrichments/1_Tables/bin/"
stagList1 = ["STAG1", "STAG2"]
stagList2 = ["69_127", "3B9_5"]

# Open File
inputFileName = fl + "1_ateb.txt"
inFile = open(inputFileName, "w")

# Stag Loop
for i in range(0,len(stagList1)):

  # Region List
  stag1 = stagList1[i]
  stag2 = stagList2[i]
  regionList = [gl+"rna.bed", gl+"tss.bed", sl+"shared.bed", sl+"STAG1_only.bed", sl+"STAG2_only.bed", cl+"CTCF.bed"]

  # Region Loop
  for region in regionList:

    # Region
    regionName = stag1 + "_" + region.split("/")[-1].split(".")[0].split("_")[0]

    # Input
    resolution = "25000"
    regionFileName = region
    abContrFileName = al + stag1 + "-.bed"
    abTreatFileName = al + stag1 + "+.bed"
    tadContrFileName = tl + stag1 + "-AUX.txt"
    tadTreatFileName = tl + stag1 + "+AUX.txt"
    #matrixContrFileName = ml + stag2 + "-.txt" # Normal
    #matrixTreatFileName = ml + stag2 + "plus.txt" # Normal
    matrixContrFileName = ml + stag2 + "-_bin.txt" # Binned
    matrixTreatFileName = ml + stag2 + "plus_bin.txt" # Binned
    temporaryLocation = tempLoc + regionName + "/"
    outputFileName = ol + regionName + ".txt"

    # Execution
    inFile.write(" ".join([resolution, regionFileName, abContrFileName, abTreatFileName, tadContrFileName, tadTreatFileName, matrixContrFileName, matrixTreatFileName, temporaryLocation, outputFileName])+"\n")

inFile.close()


