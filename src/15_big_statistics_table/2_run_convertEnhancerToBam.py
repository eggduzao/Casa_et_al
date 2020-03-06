
# Import
import os
import sys

# Input 
enhancerFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/gene_annotation/HCT116_enhancers.txt"
chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
tempLoc = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/gene_annotation/TEMP/"
outputFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/gene_annotation/HCT116_enhancers.bam"

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

# Converting enhancer to bed
enhancerFile = open(enhancerFileName, "rU")
bedFileName = tempLoc+"bedFileName.bed"
bedFile = open(bedFileName, "w")
counter = 1
for line in enhancerFile:
  ll = line.strip().split("\t")
  bedFile.write("\t".join(ll[:3]+[":".join(["E"+str(counter),ll[6],ll[8],ll[9]]), "0", "+"])+"\n")
  counter += 1
enhancerFile.close()
bedFile.close()

# Grep chromosomes
grepFileName = tempLoc+"grepFileName.bed"
command = "grep -E 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX' "+bedFileName+" > "+grepFileName
os.system(command)

# Sort bed
sortFileName = tempLoc+"sortFileName.bed"
command = "sort -k1,1 -k2,2n "+grepFileName+" > "+sortFileName
os.system(command)

# Bed To Bam
bamFileName = tempLoc+"bamFileName.bam"
command = "bedToBam -i "+sortFileName+" -g "+chromSizesFileName+" > "+bamFileName
os.system(command)

# Sort Bam
command = "samtools sort "+bamFileName+" -o "+outputFileName
os.system(command)

# Index Bam
command = "samtools index "+outputFileName
os.system(command)

# Removing all files
command = "rm -rf "+tempLoc
os.system(command)


