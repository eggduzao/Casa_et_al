
# Import
import os
import sys
from glob import glob

###################################################################################################
# INPUT
###################################################################################################

# Input
il1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/18_Process_Losada_ChIP/4_Peaks/"
su1 = "_peaks_filtered.narrowPeak"
il2 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bed_files/"
su2 = ".bed"
temporaryLocation = "/scratch/egadegu/RCT/"
outputFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/20_Comparison_SA12_Our_Losada/table.txt"

###################################################################################################
# FUNCTIONS
###################################################################################################

def file_len(fname):
  i = -1
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i + 1

def pair_compare(inputFile1Name, inputFile2Name, temporaryLocation):

  cutFileName1 = temporaryLocation + "cutFileName1.bed"
  command = "cut -f 1,2,3 "+inputFile1Name+" > "+cutFileName1
  os.system(command)

  sortFileName1 = temporaryLocation + "sortFileName1.bed"
  command = "sort -k1,1 -k2,2n "+cutFileName1+" > "+sortFileName1
  os.system(command)

  cutFileName2 = temporaryLocation + "cutFileName2.bed"
  command = "cut -f 1,2,3 "+inputFile2Name+" > "+cutFileName2
  os.system(command)

  sortFileName2 = temporaryLocation + "sortFileName2.bed"
  command = "sort -k1,1 -k2,2n "+cutFileName2+" > "+sortFileName2
  os.system(command)

  tempIntFileName = temporaryLocation + "tempIntFileName.bed"
  command = "intersectBed -u -wa -a "+sortFileName1+" -b "+sortFileName2+" > "+tempIntFileName
  os.system(command)

  tempOnly1FileName = temporaryLocation + "tempOnly1FileName.bed"
  command = "intersectBed -v -wa -a "+sortFileName1+" -b "+sortFileName2+" > "+tempOnly1FileName
  os.system(command)

  tempOnly2FileName = temporaryLocation + "tempOnly2FileName.bed"
  command = "intersectBed -v -wa -a "+sortFileName2+" -b "+sortFileName1+" > "+tempOnly2FileName
  os.system(command)

  return file_len(tempOnly1FileName), file_len(tempIntFileName), file_len(tempOnly2FileName)

###################################################################################################
# EXECUTION
###################################################################################################

# Create output location
command = "mkdir -p "+temporaryLocation
os.system(command)
outputLocation = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

# File Lists
fileList1 = ["SA1_MCF10A_ChIP-seq", "SA1_MCF10A_ChIP-seq", "SA1_MCF10A_ChIP-seq_Control", "SA1_MCF10A_ChIP-seq_Control", "SA1_MCF10A_ChIP-seq_siSA1", "SA1_MCF10A_ChIP-seq_siSA1", "SA1_MCF10A_ChIP-seq_siSA2", "SA1_MCF10A_ChIP-seq_siSA2", 
             "SA2_MCF10A_ChIP-seq", "SA2_MCF10A_ChIP-seq", "SA2_MCF10A_ChIP-seq_Control", "SA2_MCF10A_ChIP-seq_Control", "SA2_MCF10A_ChIP-seq_siSA1", "SA2_MCF10A_ChIP-seq_siSA1", "SA2_MCF10A_ChIP-seq_siSA2", "SA2_MCF10A_ChIP-seq_siSA2"]
fileList2 = ["STAG1_full_peaks", "STAG1_only", "STAG1_full_peaks", "STAG1_only", "STAG1_full_peaks", "STAG1_only", "STAG1_full_peaks", "STAG1_only", 
             "STAG2_full_peaks", "STAG2_only", "STAG2_full_peaks", "STAG2_only", "STAG2_full_peaks", "STAG2_only", "STAG2_full_peaks", "STAG2_only"]
outputFile = open(outputFileName, "w")
for i in range(0, len(fileList1)):
  fileName1 = il1 + fileList1[i] + su1
  fileName2 = il2 + fileList2[i] + su2
  a, b, c = pair_compare(fileName1, fileName2, temporaryLocation)
  outputFile.write("\t".join([fileList1[i], fileList2[i], str(a), str(b), str(c)])+"\n")
outputFile.close()


