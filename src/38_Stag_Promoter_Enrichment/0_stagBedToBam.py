
# Import
import os
import sys

# Input
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_bed_files/"
inList = ["nonpredominant", "shared", "STAG1_full_peaks", "STAG1_only", "STAG1_predominant", "STAG2_full_peaks", "STAG2_only", "STAG2_predominant"]

# Temp folder
tempLoc = "./TEMP/"
command = "mkdir -p "+tempLoc
os.system(command)

# Execution
for inName in inList:

  inBedFileName = il + inName + ".bed"
  outBamFileName = il + inName + ".bam"

  if(inName == "shared"):
    complBedFileName = tempLoc + "complBedFileName.bed"
    inBedFile = open(inBedFileName, "rU")
    complBedFile = open(complBedFileName, "w")
    counter = 1
    for line in inBedFile:
      ll = line.strip().split("\t")
      complBedFile.write("\t".join(ll+["shared"+str(counter), "1000", "+"])+"\n")
      counter += 1
    inBedFile.close()
    complBedFile.close()
  else:
    complBedFileName = tempLoc + "complBedFileName.bed"
    inBedFile = open(inBedFileName, "rU")
    complBedFile = open(complBedFileName, "w")
    for line in inBedFile:
      ll = line.strip().split("\t")
      complBedFile.write("\t".join(ll+["+"])+"\n")
    inBedFile.close()
    complBedFile.close()

  grepFileName = tempLoc + "grepFileName.bed"
  command = "grep -E 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX' "+complBedFileName+" > "+grepFileName
  os.system(command)

  sortFileName = tempLoc + "complBedFileName.bed"
  command = "sort -k1,1 -k2,2n "+grepFileName+" > "+sortFileName
  os.system(command)

  tempBamFileName = tempLoc + "tempBamFileName.bam"
  command = "bedToBam -i "+sortFileName+" -g /home/egg/rgtdata/hg19/chrom.sizes.hg19.filter > "+tempBamFileName
  os.system(command)

  command = "samtools sort "+tempBamFileName+" -o "+outBamFileName
  os.system(command)

  command = "samtools index "+outBamFileName
  os.system(command)

# Removing temp
command = "rm -rf "+tempLoc
os.system(command)


