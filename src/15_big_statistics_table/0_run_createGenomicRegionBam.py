
# Import
import os
import sys

# Input
ol = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/0_genomicRegions/"
chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
ensemblGeneFileName = "/home/egg/Projects/Wendt_Stag/Data/gene_annotation/hg19_ensembl_genes.txt"
tempLocation = "./TEMP/"
outputFileName = ol + "hg19_regions.bam"

# Execution
command = "python 0_createGenomicRegionBam.py "+" ".join([chromSizesFileName, ensemblGeneFileName, tempLocation, outputFileName])
os.system(command)


