
# Import
import os
import sys

# Input 
bedFileName = sys.argv[1]
chromSizesFileName = sys.argv[2]
tempLoc = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

# Cut
cutFileName = tempLoc+"cutFileName.bed"
command = "cut -f 1,2,3,4,5,6 "+bedFileName+" > "+cutFileName
os.system(command)

# Grep chromosomes
grepFileName = tempLoc+"grepFileName.bed"
command = "grep -E 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX' "+cutFileName+" > "+grepFileName
os.system(command)

# Put score in the name section
scoredFileName = tempLoc+"scoredFileName.bed"
grepFile = open(grepFileName, "rU")
scoredFile = open(scoredFileName, "w")
for line in grepFile:
  ll = line.strip().split()
  nameList = ll[3].split(",")
  strandList = ll[5].split(",")
  newNameList = []
  for i in range(0,len(nameList)):
    name = nameList[i].replace("::", "_")
    strand = strandList[i]
    if(name in newNameList): continue
    newNameList.append(name)
    newNameField = name+":"+ll[4]
    scoredFile.write("\t".join(ll[:3]+[newNameField,"1000",strand])+"\n")
grepFile.close()
scoredFile.close()

# Get unique entries
uniqFileName = tempLoc+"uniqFileName.bed"
command = "sort "+scoredFileName+" | uniq > "+uniqFileName
os.system(command)

# Sort bed
sortFileName = tempLoc+"sortFileName.bed"
command = "sort -k1,1 -k2,2n "+uniqFileName+" > "+sortFileName
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


