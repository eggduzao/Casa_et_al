
# Import
import os
import sys

# Gene Lists
regionFile = "/home/egg/Projects/Wendt_Stag/Results/0_Definitive_Gene_Annotation/regions.bed"
aliasFile = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
gl = "/home/egg/Projects/Wendt_Stag/Data/expression/gene_lists/"
tl = "./TEMP/"
ol = "/home/egg/Projects/Wendt_Stag/Results/11_Insulation_Plot/0_Input_Regions/"
geneList = ["downregulated_by_STAG1_only_and_specific", "downregulated_by_STAG2_only_and_specific",
            "upregulated_by_STAG1_only_and_specific", "upregulated_by_STAG2_only_and_specific"]

# Gene Loop
for geneFile in geneList:

  gg = geneFile.split("_")
  gname = gg[2] + "_" + gg[0]

  # Region List
  regionList = ["TSS", "TTS"]

  for region in regionList:

    # Input
    region = region
    aliasDictFileName = aliasFile
    regionsFileName = regionFile
    genesFileName = gl + geneFile + ".txt"
    tempLocation = tl 
    outputFileName = ol + gname + "_" + region + ".bed"

    # Execution
    command = "python 0_inputRegions.py "+" ".join([region, aliasDictFileName, regionsFileName, genesFileName, tempLocation, outputFileName])
    os.system(command)


