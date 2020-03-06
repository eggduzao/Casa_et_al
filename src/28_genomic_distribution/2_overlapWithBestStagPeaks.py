
# Import
import os
import sys
import pyBigWig
from pysam import Samfile
from numpy import isfinite, isnan

###################################################################################################
# Input
###################################################################################################

# Input
stagTotalCount = float(sys.argv[1])
ctcfTotalCount = float(sys.argv[2])
stagFileName = sys.argv[3]
stagBwSignalFileName = sys.argv[4]
ctcfSignalFileName = sys.argv[5]
genomicRegionsFileName = sys.argv[6]
chromHmmRegionsFileName = sys.argv[7]
outputFileNamePrefix = sys.argv[8]

###################################################################################################
# Functions
###################################################################################################

def fetchSignalBw(bwFile, region):
  totalSignal = 0.0
  try: valuesVec = [e if isfinite(e) and not isnan(e) else 0.0 for e in bwFile.values(region[0], region[1], region[2])]
  except Exception: return totalSignal
  return sum(valuesVec)

def fetchTotalReadsBam(bamFile, region):
  returnN = 0
  for read in bamFile.fetch(region[0], region[1], region[2]): returnN += 1
  return returnN
  
###################################################################################################
# Iterating on Stag file
###################################################################################################

# Files
stagFile = open(stagFileName, "rU")
stagBwSignalFile = pyBigWig.open(stagBwSignalFileName)
ctcfSignalFile = Samfile(ctcfSignalFileName, "rb")
genomicRegionsFile = Samfile(genomicRegionsFileName, "rb")
chromHmmRegionsFile = Samfile(chromHmmRegionsFileName, "rb")
outRegionFileName = outputFileNamePrefix+"_regions.txt"
outChromHmmFileName = outputFileNamePrefix+"_chromhmm.txt"
outRegionFile = open(outRegionFileName, "w")
outChromHmmFile = open(outChromHmmFileName, "w")

# Writing header
header = ["CHROMOSOME", "STAG_START", "STAG_END", "STAG_TOTAL_SIGNAL", "STAG_CTCF_SIGNAL", "REGION_START", "REGION_END", "REGION", "STRAND"]
outRegionFile.write("\t".join(header)+"\n")
outChromHmmFile.write("\t".join(header)+"\n")

# Iteration
for line in stagFile:

  # Initialization
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]
  
  # Fetching signal
  totalStagSignal = fetchSignalBw(stagBwSignalFile, [chrom, int(p1), int(p2)])
  totalStagSignal = str(round(totalStagSignal / ((float(p2)-float(p1))/1000 * stagTotalCount/1000000),4))
 
  # Fetching CTCF signal
  totalCtcfSignal = fetchTotalReadsBam(ctcfSignalFile, [chrom, int(p1), int(p2)])
  totalCtcfSignal = str(round(totalCtcfSignal / ((float(p2)-float(p1))/1000 * ctcfTotalCount/1000000),4))
  
  # Fetching chromhmm states
  chromHmmFetch = chromHmmRegionsFile.fetch(chrom, int(p1), int(p2))
  for read in chromHmmFetch:
    regionStart = read.pos
    regionEnd = read.aend
    region = read.qname
    outChromHmmFile.write("\t".join([chrom, p1, p2, totalStagSignal, str(regionStart), str(regionEnd), region, "."])+"\n")

  # Fetching regions
  regions = []
  regionFetch = genomicRegionsFile.fetch(chrom, int(p1), int(p2))
  for read in regionFetch:
    regionStart = read.pos
    regionEnd = read.aend
    region = read.qname
    rr = region.split(":")
    if(rr[0] == "TTS" or rr[0] == "INTRON" or rr[0] == "EXON"): continue
    if(rr[0] == "ENHANCER" or rr[0] == "SUPERENHANCER" or rr[0] == "INTERGENIC" or rr[0] == "CTCF"): strand = "."
    else:
      strand = "+"
      if(read.is_reverse): strand = "-"
    outRegionFile.write("\t".join([chrom, p1, p2, totalStagSignal, totalCtcfSignal, str(regionStart), str(regionEnd), region, strand])+"\n")

# Closing files
stagFile.close()
stagBwSignalFile.close()
ctcfSignalFile.close()
genomicRegionsFile.close()
chromHmmRegionsFile.close()
outRegionFile.close()
outChromHmmFile.close()


