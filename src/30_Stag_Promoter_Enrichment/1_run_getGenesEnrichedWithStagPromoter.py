
# Import
import os
import sys

# Peak List
aliasFile = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
regionFile = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/1_genomic_regions/regions.bed"
el = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/expression/"
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_bed_files/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/1_stag_gene_list/"
peakList = ["STAG1bestpeaks_filter", "STAG2bestpeaks_filter", "STAG1only_mrg_replicates", "STAG2only_mtg_replicates"]

# Peak Loop
for peakName in peakList:

  # Parameters
  name = peakName.split("_")[0]
  stag = "STAG2"
  if("STAG1" in name): stag = "STAG1"

  # Input
  peakExt = "1000"
  aliasFileName = aliasFile
  expressionFileName = el + stag + "_minus.txt"
  regionFileName = regionFile
  stagPeakFileName = il + peakName + ".bam"
  outputActiveFileName = ol + name + "_active.txt"
  outputInactiveFileName = ol + name + "_inactive.txt"

  # Execution
  command = "python 1_getGenesEnrichedWithStagPromoter.py "+" ".join([peakExt, aliasFileName, expressionFileName, regionFileName, stagPeakFileName, outputActiveFileName, outputInactiveFileName])
  os.system(command)


