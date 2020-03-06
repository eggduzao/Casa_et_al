
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
aliasFileName = sys.argv[1]
geneFileName = sys.argv[2]
treatFileName = sys.argv[3]
controlFileName = sys.argv[4]
tempLoc = sys.argv[5]
outputFileName = sys.argv[6]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def createVectorOfGenes(expFileName, aliasDict, geneDict):
  geneExpVec = []
  expFile = open(expFileName, "rU")
  for line in expFile:
    ll = line.strip().split("\t")
    g = ll[0]; exp = float(ll[1])
    try: gg = aliasDict[g.upper()]
    except Exception: gg = g.upper()
    try: gene = geneDict[gg]
    except Exception: continue
    geneLen = float(int(gene[2]) - int(gene[1]))
    exp = str(exp / geneLen)
    geneExpVec.append([gg, exp])
  expFile.close()
  return geneExpVec

###################################################################################################
# Gene Dictionary
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

# Reading genes
geneDict = dict()
geneFile = open(geneFileName, "rU")
for line in geneFile:
  ll = line.strip().split("\t")
  try: gene = aliasDict[ll[3].upper()]
  except Exception: gene = ll[3].upper()
  geneDict[gene] = ll
geneFile.close()
  
###################################################################################################
# Creating table
###################################################################################################

# Treatment expression
treatVec = createVectorOfGenes(treatFileName, aliasDict, geneDict)

# Control expression
controlVec = createVectorOfGenes(controlFileName, aliasDict, geneDict)

# Writing vectors
outputFile = open(outputFileName, "w")
outputFile.write("\t".join(["GENE", "EXP", "TYPE"])+"\n")
for vec in treatVec: outputFile.write("\t".join(vec+["T"])+"\n")
for vec in controlVec: outputFile.write("\t".join(vec+["C"])+"\n")
outputFile.close()

# Termination
command = "rm -rf "+tempLoc
os.system(command)


