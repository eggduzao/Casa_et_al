
# Import
import os
import sys

# Peak List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/35_Venn_SA12_Rad21_Ctcf_Smc3/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/35_Venn_SA12_Rad21_Ctcf_Smc3/0_Peaks/"
tl = "/scratch/egadegu/COT/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/35_Venn_SA12_Rad21_Ctcf_Smc3/1_Tables/"
peakList = ["STAG1", "STAG2", "SMC3", "CTCF", "RAD21"]

# Open File
inputFileName = fl + "1_cot.txt"
inFile = open(inputFileName, "w")

# Input
halfExt = "100"
peakFileList = ",".join([il + e + ".bed" for e in peakList])
temporaryLocation = tl
outputPairFileName = ol + "table_pair.txt"
outputTripleFileName = ol + "table_triple.txt"

# Execution
inFile.write(" ".join([halfExt, peakFileList, temporaryLocation, outputPairFileName, outputTripleFileName])+"\n")

# Peak List
tl = "/scratch/egadegu/COT2/"
peakList = ["STAG1_REP1", "STAG1_REP2", "STAG2_REP1", "STAG2_REP2"]

# Input
halfExt = "100"
peakFileList = ",".join([il + e + ".bed" for e in peakList])
temporaryLocation = tl
outputPairFileName = ol + "table_pair_rep.txt"
outputTripleFileName = ol + "table_triple_rep.txt"

# Execution
inFile.write(" ".join([halfExt, peakFileList, temporaryLocation, outputPairFileName, outputTripleFileName])+"\n")

inFile.close()


