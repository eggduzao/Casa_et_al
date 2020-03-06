
# Import
import os
import sys

# Expression List
counter = 1
aliasFile = "/home/egusmao/rgtdata/hg19/alias_human_booster.txt"
regionFile = "/projects/ag-papan/eduardo/Papantonis_Stag/Previous_Results/28_genomic_distribution/1_genomic_regions/regions.bed"
inf = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/7_Stag_Promoter_Enrichment/input/"
el = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/expression/filtered_all/"
rl = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_bed_files/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/7_Stag_Promoter_Enrichment/1_Promoter_Enrichment_Table/"
expressionList = ["STAG1_all", "STAG2_all"]

# Expression Loop
for exp in expressionList:

  # Region List
  regionList = ["STAG1_only", "STAG2_only", "shared", "STAG1_predominant", "STAG2_predominant", "nonpredominant"]

  # Region Loop
  for reg in regionList:

    # Parameters
    region = reg
    expression = exp.split("_")[0]

    # Input
    randomRepeats = "1000"
    promoterExt = "0"
    aliasFileName = aliasFile
    genomicRegionsFileName = regionFile
    stagRegionsFileName = rl + reg + ".bam"
    expressionFileName = el + exp + ".txt"
    outputTableFileName = ol + "exp_" + expression + "_region_" + region + ".txt"

    # Execution
    inFile = open(inf+str(counter)+"_pet.txt", "w")
    inFile.write("\n".join([randomRepeats, promoterExt, aliasFileName, genomicRegionsFileName, stagRegionsFileName, expressionFileName, outputTableFileName]))
    inFile.close()
    counter += 1


