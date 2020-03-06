
# Import
import os
import sys

# Loop List
il1 = "/home/egg/Projects/Wendt_Stag/Previous_Results/10_All_Differential_Loopograms/input/10K_norm/"
il2 = "/home/egg/Projects/Wendt_Stag/Previous_Results/3_Loopograms/matrix_new/10K_norm/"
ol = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/"
loopList = ["R_STAG1_RM_STAG1_intersection", "R_STAG1_RM_STAG1_minusOnly", "R_STAG1_RM_STAG1_plusOnly", "69_127plus.txt:STAG1_regions", "69_127-.txt:STAG1_regions"]

# Loop Loop
for loopName in loopList:

  # Parameters
  il = il1
  suff = ".bed"
  if("txt" in loopName):
    il = il2
    suff = ".txt"

  # Motif List
  ml = "/home/egg/Projects/Wendt_Stag/Data/tf/motifs/"
  motifList = ["HCT116_ChIP-seq_CTCF_HAIB", "HCT116_ChIP-seq_CTCF_UW", "HCT116_ChIP-seq_CTCF_BROAD"]

  # Motif Loop
  for motifName in motifList:

    # Parameters
    name = "__".join([loopName, "_".join(motifName.split("_")[2:])])

    # Input
    loopFileName = il + loopName + suff
    motifFileName = ml + motifName + ".bed"
    tempLocation = "./TEMP/" + name + "/"
    outputFileName = ol + name + ".txt"

    # Execution
    command = "python 1_ctcfOrientation.py "+" ".join([loopFileName, motifFileName, tempLocation, outputFileName])
    os.system(command)


