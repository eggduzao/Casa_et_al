
# Import
import os
import sys
from pysam import Samfile

###################################################################################################
# Input
###################################################################################################

# Input
motifFileNameList = sys.argv[1].split(",")
chromSizesFileName = sys.argv[2]
summitFileName = sys.argv[3]
tempLocation = sys.argv[4]
outputFileName = sys.argv[5]

# Initialization
summitExt = 50
command = "mkdir -p "+tempLocation
os.system(command)

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

###################################################################################################
# Functions
###################################################################################################

# Convert bed to bam
def bedToBam(bedFileName, chromSizesFileName, finalBamFileName, tempLoc):

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
  command = "samtools sort "+bamFileName+" -o "+finalBamFileName
  os.system(command)

  # Index Bam
  command = "samtools index "+finalBamFileName
  os.system(command)

  # Return
  return 0

###################################################################################################
# Execution
###################################################################################################

# Convert motif files to bam
newMotifFileNameList = []
for motifFileName in motifFileNameList:
  newBamFile = ".".join(motifFileName.split(".")[:-1])+".bam"
  bedToBam(motifFileName, chromSizesFileName, newBamFile, tempLocation)
  newMotifFileNameList.append(newBamFile)

# Opening all motif files
motifFileList = [Samfile(e, "rb") for e in newMotifFileNameList]

# Iterating on summit file
summitFile = open(summitFileName, "rU")
outputMotifsFileName = tempLocation+"outputMotifsFileName.bed"
outputMotifsFile = open(outputMotifsFileName, "w")
for line in summitFile:

  # Initialization
  ll = line.strip().split("\t")
  if(ll[0] not in chrList): continue

  # Coordinates
  chrom = ll[0]
  mid = (int(ll[1])+int(ll[2]))/2
  p1 = mid - summitExt; p2 = mid + summitExt

  # Fetching motifs
  bestMotif = []
  bestScore = -999999.
  for motifFile in motifFileList:
    try: motifFetch = motifFile.fetch(chrom, p1, p2)
    except Exception: continue
    for read in motifFetch:
      rr = read.qname.split(":")
      motifName = rr[0]
      motifScore = rr[1]
      if(float(motifScore) <= bestScore): continue
      bestScore = float(motifScore)
      startx = str(read.pos)
      endx = str(read.aend)
      strand = "+"
      if(read.is_reverse): strand = "-"
      bestMotif = [chrom, startx, endx, motifName, motifScore, strand]
  
  # Writing best motif
  if(len(bestMotif) > 1): outputMotifsFile.write("\t".join(bestMotif)+"\n")

# Closing files
summitFile.close()
outputMotifsFile.close()
for e in motifFileList: e.close()

# Sorting output file
command = "sort -k1,1 -k2,2n "+outputMotifsFileName+" > "+outputFileName
os.system(command)

# Convert outputFileName to bam
outputFileBamName = ".".join(outputFileName.split(".")[:-1])+".bam"
bedToBam(outputFileName, chromSizesFileName, outputFileBamName, tempLocation)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


"""

# Import
import os
import sys

###################################################################################################
# Input
###################################################################################################

# Input
motifFileNameList = sys.argv[1].split(",")
summitFileName = sys.argv[2]
tempLocation = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
summitExt = 100
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Function
###################################################################################################

def fileLen(fileName):
  if(os.stat(fileName).st_size <= 0): return 0
  i = 0
  with open(fileName) as f:
    for i, l in enumerate(f): pass
  return i + 1

###################################################################################################
# Execution
###################################################################################################

# Fetching motifs
if(len(motifFileNameList) > 1):
  unsortedMainMotifFile = tempLocation+"unsortedMainMotifFile.bed"
  command = "cat "+motifFileNameList[0]+" > "+unsortedMainMotifFile
  os.system(command)  
  for mfn in motifFileNameList[1:]: 
    newNIntFileName = tempLocation+mfn.split("/")[-1].split(".")[0]+"_newNIntFileName.bed"
    command = "intersectBed -wa -v -a "+mfn+" -b "+motifFileNameList[0]+" > "+newNIntFileName
    os.system(command)  
    command = "cat "+newNIntFileName+" >> "+unsortedMainMotifFile
    os.system(command) 
  sortedMainMotifFile = tempLocation+"sortedMainMotifFile.bed"
  command = "sort -k1,1 -k2,2n "+unsortedMainMotifFile+" > "+sortedMainMotifFile
  os.system(command)
  motifFileName = tempLocation+"motifFileName.bed"
  command = "mergeBed -c 4,5,6 -o collapse,mean,collapse -i "+sortedMainMotifFile+" > "+motifFileName
  os.system(command)
else: motifFileName = motifFileNameList[0]

# Creating peak file
summitFile = open(summitFileName, "rU")
peakFileName = tempLocation+"peakFileName.bed"
peakFile = open(peakFileName, "w")
for line in summitFile:
  ll = line.strip().split("\t")
  mid = (int(ll[1])+int(ll[2]))/2
  peakFile.write("\t".join([ll[0],str(mid-summitExt),str(mid+summitExt)])+"\n")
summitFile.close()
peakFile.close()

# Reducing motifFile
reducedMotifFileName = tempLocation+"reducedMotifFileName.bed"
command = "intersectBed -wa -u -a "+motifFileName+" -b "+peakFileName+" > "+reducedMotifFileName
os.system(command)

# Fetching only the summit of the peaks and higest scoring motif it overlaps
motifList = []
peakFile = open(peakFileName, "rU")
for line in peakFile:
  tempPeakFileName = tempLocation+"tempPeakFileName.bed"
  tempPeakFile = open(tempPeakFileName, "w")
  tempPeakFile.write(line)
  tempPeakFile.close()
  intersectFileName = tempLocation+"intersectFileName.bed"
  command = "intersectBed -wa -u -a "+reducedMotifFileName+" -b "+tempPeakFileName+" > "+intersectFileName
  os.system(command)
  if(fileLen(intersectFileName) <= 0): continue
  intersectFile = open(intersectFileName, "rU")
  currMotif = None
  currScore = -999999.9
  for gline in intersectFile:
    gg = gline.strip().split("\t")
    if(float(gg[4]) > currScore):
      currScore = float(gg[4])
      currMotif = gg
  intersectFile.close()
  if(currMotif): motifList.append(currMotif)
peakFile.close()

# Writing motifs
outputFile = open(outputFileName, "w")
for vec in motifList: outputFile.write("\t".join(vec)+"\n")
outputFile.close()

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)
"""


