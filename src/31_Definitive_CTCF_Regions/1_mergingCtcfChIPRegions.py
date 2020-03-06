
# Import
import os
import sys
from pysam import Samfile

# Input
chromSizesFileName = sys.argv[1]
chipRegionFileName = sys.argv[2]
motifBamFileName = sys.argv[3]
tempLoc = sys.argv[4]
outputBedFileName = sys.argv[5]
outputBamFileName = sys.argv[6]

# Creating tempLoc
chromList = ["chr"+str(e) for e in range(1,23)+["X"]]
command = "mkdir -p "+tempLoc
os.system(command)

# Iterating throught the regions
chipRegionFile = open(chipRegionFileName, "rU")
motifBamFile = Samfile(motifBamFileName, "rb")
tempBedFileName = tempLoc + "tempBedFileName"
tempBedFile = open(tempBedFileName, "w")
for line in chipRegionFile:
  ll = line.strip().split("\t")
  chrom = ll[0]; start = int(ll[1]); end = int(ll[2])
  if(chrom not in chromList): continue
  mfetch = motifBamFile.fetch(chrom, start, end)
  bestScore = -9999
  bestMotif = None
  for read in mfetch:
    reference_start = str(read.reference_start)
    reference_end = str(read.reference_end)
    rr = read.query_name.split(":")
    name = rr[0]
    score = float(rr[1])
    strand = "+"
    if(read.is_reverse): strand = "-"
    if(score > bestScore):
      bestScore = score
      bestMotif = [chrom, reference_start, reference_end, name, str(score), strand]
  if(bestMotif): tempBedFile.write("\t".join(bestMotif)+"\n")
chipRegionFile.close()
motifBamFile.close()
tempBedFile.close()

# Sort motifs
command = "sort -k1,1 -k2,2n "+tempBedFileName+" > "+outputBedFileName
os.system(command)

# Create bam
tempBamFileName = tempLoc + "tempBamFileName.bed"
command = "bedToBam -i "+outputBedFileName+" -g "+chromSizesFileName+" > "+tempBamFileName
os.system(command)

# Sort bam
command = "samtools sort "+tempBamFileName+" -o "+outputBamFileName
os.system(command)

# Index bam
command = "samtools index "+outputBamFileName
os.system(command)

# Removing temp loc
command = "rm -rf "+tempLoc
os.system(command)


