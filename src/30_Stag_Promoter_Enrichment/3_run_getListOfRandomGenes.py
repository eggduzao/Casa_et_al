
# Import
import os
import sys

# Peak List
gl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/2_active_inactive_genes/"
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/1_stag_gene_list/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/3_control_gene_list/"
peakList = ["STAG1bestpeaks_filter", "STAG2bestpeaks_filter", "STAG1only_mrg_replicates", "STAG2only_mtg_replicates"]

# Peak Loop
for peakName in peakList:

  # Parameters
  name = peakName.split("_")[0]
  stag = "STAG1"
  if("STAG2" in name): stag = "STAG2"

  # Input
  inputActiveFileName = il + name + "_active.txt"
  inputInactiveFileName = il + name + "_inactive.txt"
  activeGenesFileName = gl + stag + "_active.txt"
  inactiveGenesFileName = gl + stag + "_inactive.txt"
  tempLocation = ol + "TEMP/"
  outputActiveControlFileName = ol + name + "_active.txt"
  outputInactiveControlFileName = ol + name + "_inactive.txt"

  # Execution
  command = "python 3_getListOfRandomGenes.py "+" ".join([inputActiveFileName, inputInactiveFileName, activeGenesFileName, inactiveGenesFileName, tempLocation, outputActiveControlFileName, outputInactiveControlFileName])
  os.system(command)


