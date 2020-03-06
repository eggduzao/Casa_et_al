
# Import
import os
import sys

# Input 
chromFileName = "/home/egg/Projects/Papantonis_Stag/Results/16_chromhmm/HCT116_20_dense.bed"
chromSizesFileName = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
tempLoc = "/home/egg/Projects/Papantonis_Stag/Results/16_chromhmm/TEMP/"
outputFileName = "/home/egg/Projects/Papantonis_Stag/Results/16_chromhmm/HCT116_20_dense.bam"

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

# Annotation dictionary
annDict = {"1": "CTCF/INSULATOR",
           "2": "CTCF/INSULATOR",
           "3": "POISED_ENHANCER",
           "4": "WEAK_ENHANCER",
           "5": "STRONG_ENHANCER",
           "6": "STRONG_ENHANCER",
           "7": "WEAK_ENHANCER",
           "8": "INACTIVE_GENE_BODY",
           "9": "ELONGATION",
           "10": "INACTIVE_GENE_BODY",
           "11": "ELONGATION",
           "12": "ACTIVE_TSS",
           "13": "POISED_ENHANCER",
           "14": "WEAK_PROMOTER",
           "15": "STRONG_PROMOTER",
           "16": "CONSTITUTIVE_HETEROCHROMATIN",
           "17": "OTHER",
           "18": "REPRESSED_POLYCOMB",
           "19": "CONSTITUTIVE_HETEROCHROMATIN",
           "20": "CONSTITUTIVE_HETEROCHROMATIN"}

# Annotating chromhmm
chromFile = open(chromFileName, "rU")
bedFileName = tempLoc+"bedFileName.bed"
bedFile = open(bedFileName, "w")
chromFile.readline()
for line in chromFile:
  ll = line.strip().split("\t")
  bedFile.write("\t".join(ll[:3]+[annDict[ll[3]], "0", "+"])+"\n")
chromFile.close()
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


