
# Import
import os
import sys

# Stag List
aliasFile = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
geneFile = "/home/egg/rgtdata/hg19/genes_RefSeq_hg19.bed"
ilt = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/1_stag_gene_list/"
ilc = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/3_control_gene_list/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/4_expression_table/"
stagList = ["STAG1", "STAG2"]

# Stag Loop
for stag in stagList:

  # Input
  aliasFileName = aliasFile
  geneFileName = geneFile
  treatFileName = ilt + stag + "only_active.txt"
  controlFileName = ilc + stag + "only_active.txt"
  tempLoc = ol + "TEMP/"
  outputFileName = ol + stag + "only_active.txt"

  # Execution
  command = "python 4_expressionTable.py "+" ".join([aliasFileName, geneFileName, treatFileName, controlFileName, tempLoc, outputFileName])
  os.system(command)


