
# Import
import os
import sys

# Input
aliasFileName = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
ensemblDictFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/gene_annotation/hg19_ensembl_dictionary_filtered.txt"
ensemblGeneFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/gene_annotation/hg19_ensembl_genes.txt"
ensemblBamFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/gene_annotation/hg19_ensembl_genes.bam"
ensemblNoIEBedFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/gene_annotation/hg19_ensembl_genes_noExonIntron.bed"
ensemblNoIEBamFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/gene_annotation/hg19_ensembl_genes_noExonIntron.bam"
enhancerFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/gene_annotation/HCT116_enhancers_filter.bam"
h3k4me1FileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/histones/macs/HCT116_ChIP-seq_H3K4me1_USC_peaks_filter.bam"
h3k4me3FileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/histones/macs/HCT116_ChIP-seq_H3K4me3_UW_peaks_filter.bam"
h3k27acFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/histones/macs/HCT116_ChIP-seq_H3K27ac_USC_peaks_filter.bam"
ctcfFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/tf/macs/HCT116_ChIP-seq_CTCF_BROAD_peaks_filter.bam"
tempLocation = "/home/egg/Projects/Wendt_Stag/Code/0_Definitive_Gene_Annotation/TEMP/"
outputBedFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/0_Definitive_Gene_Annotation/regions.bed"
outputBamFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/0_Definitive_Gene_Annotation/regions.bam"

# Execution
command = "python 1_createGenomicRegionFiles.py "+" ".join([aliasFileName, chromSizesFileName, ensemblDictFileName, ensemblGeneFileName, ensemblBamFileName, ensemblNoIEBedFileName, ensemblNoIEBamFileName, enhancerFileName, h3k4me1FileName, h3k4me3FileName, h3k27acFileName, ctcfFileName, tempLocation, outputBedFileName, outputBamFileName])
os.system(command)


