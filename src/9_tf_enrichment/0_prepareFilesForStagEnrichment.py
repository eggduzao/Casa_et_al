
# Import
import os
import sys

# Input
ext = int(sys.argv[1])
peakFileName = sys.argv[2]
stagFileName = sys.argv[3]
tempLocation = sys.argv[4]
outputFileName = sys.argv[5]

# Initialization
inName = peakFileName.split("/")[-1].split(".")[0]
inExt = peakFileName.split("/")[-1].split(".")[-1]
command = "mkdir -p "+tempLocation
os.system(command)

# Preparing peaks
if(inExt == "narrowPeak"):
  newPeakFileName = tempLocation+inName+"_ext.bed"
  inFile = open(peakFileName,"rU")
  extFile = open(newPeakFileName,"w")
  for line in inFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); summit = int(float(ll[9]))
    mid = p1 + summit
    e1 = str(mid - ext); e2 = str(mid + ext)
    extFile.write("\t".join([chrom, e1, e2])+"\n")
  inFile.close()
  extFile.close()

# Preparing summits
elif("summits" in inName):
  newPeakFileName = tempLocation+inName+"_ext.bed"
  inFile = open(peakFileName,"rU")
  extFile = open(newPeakFileName,"w")
  for line in inFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])
    mid = (int(p1) + int(p2)) / 2
    e1 = str(mid - ext); e2 = str(mid + ext)
    extFile.write("\t".join([chrom, e1, e2])+"\n")
  inFile.close()
  extFile.close()

# Preparing footprints
else:
  newPeakFileName = tempLocation+inName+"_ext.bed"
  inFile = open(peakFileName,"rU")
  extFile = open(newPeakFileName,"w")
  for line in inFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])
    mid = (int(p1) + int(p2)) / 2
    e1 = str(mid - ext); e2 = str(mid + ext)
    extFile.write("\t".join([chrom, e1, e2])+"\n")
  inFile.close()
  extFile.close()

# Intersection with STAG file
command = "intersectBed -wa -u -a "+newPeakFileName+" -b "+stagFileName+" > "+outputFileName
os.system(command)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


