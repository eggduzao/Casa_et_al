
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

# Input
peakExt = int(sys.argv[1])
aliasFileName = sys.argv[2]
expressionFileName = sys.argv[3]
regionFileName = sys.argv[4]
stagPeakFileName = sys.argv[5]
outputActiveFileName = sys.argv[6]
outputInactiveFileName = sys.argv[7]

###################################################################################################
# Functions
###################################################################################################

def check_bam_at_least_one_read(bam_file, region):
  res = False
  for read in bam_file.fetch(region[0], region[1], region[2]):
    res = True
    break
  return res
  
###################################################################################################
# Expression Dictionary
###################################################################################################

# Reading alias
aliasDict = dict() # alias -> gene_symbol
aliasFile = open(aliasFileName,"rU")
for line in aliasFile:
  ll = line.strip().split("\t")
  value = ll[1]
  geneList = [ll[0],ll[1]]+ll[2].split("&")
  for g in geneList: aliasDict[g.upper()] = value.upper()
aliasFile.close()

# Expression dictionary
expDict = dict()
expressionFile = open(expressionFileName, "rU")
for line in expressionFile:
  ll = line.strip().split("\t")
  g = ll[0].upper()
  value = ll[1]
  try: gene = aliasDict[g]
  except Exception: gene = g
  expDict[gene] = value
expressionFile.close()

###################################################################################################
# Fetching STAG Genes
###################################################################################################

# Opening files
regionFile = open(regionFileName, "rU")
stagPeakFile = Samfile(stagPeakFileName, "rb")
outputActiveFile = open(outputActiveFileName, "w")
outputInactiveFile = open(outputInactiveFileName, "w")

# Iterating on region file
for line in regionFile:

  ll = line.strip().split("\t")
  
  chromosome = ll[0]; start = ll[1]; end = ll[2]; regionList = ll[3].split(":"); score = ll[4]; strand = ll[5]
  if(regionList[0] != "PROMOTER"): continue

  check = check_bam_at_least_one_read(stagPeakFile, [chromosome, int(start)-peakExt, int(end)+peakExt])
  if(not check): continue

  try: gene = aliasDict[regionList[1].upper()]
  except Exception: gene = regionList[1].upper()
  try: exp = expDict[gene]
  except Exception: continue

  if(regionList[2] == "ACTIVE"): outputActiveFile.write("\t".join([gene, exp])+"\n")
  elif(regionList[2] == "INACTIVE"): outputInactiveFile.write("\t".join([gene, exp])+"\n")

# Termination
regionFile.close()
stagPeakFile.close()
outputActiveFile.close()
outputInactiveFile.close()


