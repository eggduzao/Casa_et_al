
# Import
import os
import sys

# Input
chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
chipRegionFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/31_Definitive_CTCF_Regions/0_input_files/CTCF_merged_chip.bed"
motifBamFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/31_Definitive_CTCF_Regions/0_input_files/CTCF_merged_motifs.bam"
tempLoc = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/31_Definitive_CTCF_Regions/1_CtcfChIPRegions/TEMP/"
outputBedFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/31_Definitive_CTCF_Regions/1_CtcfChIPRegions/CTCF_regions.bed"
outputBamFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/31_Definitive_CTCF_Regions/1_CtcfChIPRegions/CTCF_regions.bam"

# Execution
command = "python 1_mergingCtcfChIPRegions.py "+" ".join([chromSizesFileName, chipRegionFileName, motifBamFileName, tempLoc, outputBedFileName, outputBamFileName])
os.system(command)


