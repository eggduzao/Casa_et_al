
# Import
import os
import sys

###################################################################################################
# Fetching all genes
###################################################################################################

# Input
aliasFileName = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
regionsFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Previous_Results/28_genomic_distribution/1_genomic_regions/regions.bed"

# Creating alias dictionary
aliasDict = dict()
aliasFile = open(aliasFileName,"rU")
for line in aliasFile:
  ll = line.strip().split("\t")
  value = ll[1]
  geneList = [ll[0],ll[1]]+ll[2].split("&")
  for g in geneList: aliasDict[g.upper()] = value.upper()
aliasFile.close()

# Creating dictionary with all genes
allGenesDict = dict() # GENE_NAME => [CHROMOSOME, START, END, ACTIVE_STATUS, GENE_LENGTH, STRAND]
regionsFile = open(regionsFileName, "rU")
for line in regionsFile:
  ll = line.strip().split("\t")
  chrom = ll[0]; start = ll[1]; end = ll[2]; nn = ll[3].split(":"); length = str(int(end)-int(start)); strand = ll[5]
  region = nn[0]; name = nn[1]; activeStatus = nn[2]
  if(region != "GENE"): continue
  gene = name
  try: gene = aliasDict[name.upper()]
  except Exception: gene = name.upper()
  toWrite = [chrom, start, end, activeStatus, length, strand]
  allGenesDict[gene] = toWrite
regionsFile.close()

###################################################################################################
# Fetching final region + expression files
###################################################################################################

# Parameters
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/expression/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/8_Meta_Plots/0_Regions/"
statusList = ["filtered_all", "filtered_up", "filtered_down", "filtered_miss"]
stagList = ["STAG1_all", "STAG2_all"]

# Status Loop
for status in statusList:

  # Stag Loop
  for stag in stagList:

    # Input
    exprFileName = il + status + "/" + stag + ".txt"
    outputFileName = ol + stag.split("_")[0] + "_" + status.split("_")[1] + ".bed"
    # BED = CHROM, START(Gene), END(Gene), GENE_NAME, LOG2_EXP_FC, STRAND.

    # Creating final file
    exprFile = open(exprFileName, "rU")
    outputFile = open(outputFileName, "w")
    exprFile.readline()
    for line in exprFile:
      ll = line.strip().split("\t")
      name = ll[1]; log2fc = ll[2]
      try: gene = aliasDict[name.upper()]
      except Exception: gene = name.upper()
      try: geneVec = allGenesDict[gene]
      except Exception:
        print outputFileName, gene
        continue
      chrom = geneVec[0]; start = geneVec[1]; end = geneVec[2]; strand = geneVec[5]
      toWrite = [chrom, start, end, gene, log2fc, strand]
      outputFile.write("\t".join(toWrite)+"\n")
    exprFile.close()
    outputFile.close()


