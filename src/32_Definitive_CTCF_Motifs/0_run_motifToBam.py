
# Import
import os
import sys

# Input
chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
inFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/31_Definitive_CTCF_Regions/0_input_files/CTCF_merged_motifs.bed"
outFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/31_Definitive_CTCF_Regions/0_input_files/CTCF_merged_motifs.bam"
tempLoc = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/31_Definitive_CTCF_Regions/0_input_files/TEMP/"

# Creating tempLoc
command = "mkdir -p "+tempLoc
os.system(command)

# Correct motif file
correctFileName = tempLoc + "correctFileName.bed"
inFileName = open(inFileName, "rU")
correctFile = open(correctFileName, "w")
for line in inFileName:
  ll = line.strip().split("\t")
  correctFile.write("\t".join(ll[:3]+[":".join(ll[3:5]), "1000", ll[5]])+"\n")
inFileName.close()
correctFile.close()

# Grep chromosomes
grepFileName = tempLoc + "grepFileName.bed"
command = "grep -E 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX' "+correctFileName+" > "+grepFileName
os.system(command)

# Sort motifs
sortFileName = tempLoc + "sortFileName.bed"
command = "sort -k1,1 -k2,2n "+grepFileName+" > "+sortFileName
os.system(command)

# Create bam
tempBamFileName = tempLoc + "tempBamFileName.bed"
command = "bedToBam -i "+sortFileName+" -g "+chromSizesFileName+" > "+tempBamFileName
os.system(command)

# Sort bam
command = "samtools sort "+tempBamFileName+" -o "+outFileName
os.system(command)

# Index bam
command = "samtools index "+outFileName
os.system(command)

# Removing temp loc
command = "rm -rf "+tempLoc
os.system(command)


