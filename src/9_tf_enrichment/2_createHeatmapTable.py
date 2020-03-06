
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
aliasFileName = sys.argv[1]
expressionListFileName1 = sys.argv[2]
expressionListFileName2 = sys.argv[3]
inputMeTableName1 = sys.argv[4]
inputMeTableName2 = sys.argv[5]
outputFileName = sys.argv[6]

###################################################################################################
# Functions
###################################################################################################

# Table to dictionary
def tableToDict(inputFileName, expDict, aliasDict):
  resDict = dict()
  inputFile = open(inputFileName, "rU")
  inputFile.readline()
  for line in inputFile:
    ll = line.strip().split("\t")
    if("_HUMAN" in ll[0]): 
      factor = ll[0].split("_")[0]
      if(aliasDict):
        try:
          factor = aliasDict[factor]
          im_so_handsome = expDict[factor]
        except Exception: continue
    else:
      factor = ll[0].split("(")[0].split(".")[-1].upper()
      factorList = factor.split("::")
      if(aliasDict):
        factorListNew = []
        for f in factorList:
          try:
            f2 = aliasDict[factor]
            im_so_handsome = expDict[f2]
          except Exception: continue
          factorListNew.append(f2)
      else: factorListNew = factorList
      factor = "c".join(factorList)
    try:
      if(float(resDict[factor]) > float(ll[1])): resDict[factor] = ll[1]
    except Exception: resDict[factor] = ll[1]
  inputFile.close()
  return resDict

# Read expression dictionary
def readExpDict(expressionListFileName):
  expressionListFile = open(expressionListFileName, "rU")
  expDict = dict()
  for line in expressionListFile: expDict[line.strip()] = True
  expressionListFile.close()
  return expDict

###################################################################################################
# Execution
###################################################################################################

# Alias dictionary
aliasDict = None
if(aliasFileName != "."):
  aliasDict = dict()
  aliasFile = open(aliasFileName,"rU")
  for line in aliasFile:
    ll = line.strip().split("\t")
    value = ll[1]
    geneList = [ll[0],ll[1]]+ll[2].split("&")
    for g in geneList: aliasDict[g] = value
  aliasFile.close()

# Expression dict
expDict1 = None
expDict2 = None
if(expressionListFileName1 != "."): expDict1 = readExpDict(expressionListFileName1)
if(expressionListFileName2 != "."): expDict2 = readExpDict(expressionListFileName2)

# Reading dictionaries
dict1 = tableToDict(inputMeTableName1, expDict1, aliasDict)
dict2 = tableToDict(inputMeTableName2, expDict2, aliasDict)

# Fetching all dictionary keys
keys1 = dict1.keys()
keys2 = dict2.keys()
allKeys = sorted(list(set(keys1).union(set(keys2))))

# Writing output
outputFile = open(outputFileName, "w")
header = ["FACTOR", "COND1", "COND2"]
outputFile.write("\t".join(header)+"\n")
for k in allKeys:
  try: value1 = dict1[k]
  except Exception: value1 = "1.0"
  try: value2 = dict2[k]
  except Exception: value2 = "1.0"
  outputFile.write("\t".join([k, value1, value2])+"\n")
outputFile.close()


