
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
sys.path = ["/home/egusmao/.local/lib/python2.7/site-packages"] + sys.path
from pysam import Samfile
from random import randint, seed

# Input
randomRepeats = int(sys.argv[1])
promoterExt = int(sys.argv[2])
aliasFileName = sys.argv[3]
genomicRegionsFileName = sys.argv[4]
stagRegionsFileName = sys.argv[5]
expressionFileName = sys.argv[6]
outputTableFileName = sys.argv[7]

# Initialization
seed(111)

###################################################################################################
# Functions
###################################################################################################

def check_bam_at_least_one_read(bam_file, region):
  res = False
  for read in bam_file.fetch(region[0], region[1], region[2]):
    res = True
    break
  return res

def get_alias_dictionary(alias_file_name):

  # Creating alias dictionary
  aliasDict = dict() # ALIAS -> GENE_SYMBOL
  aliasFile = open(alias_file_name,"rU")
  for line in aliasFile:
    ll = line.strip().split("\t")
    value = ll[1]
    geneList = [ll[0],ll[1]]+ll[2].split("&")
    for g in geneList: aliasDict[g.upper()] = value.upper()
  aliasFile.close()

  # Returning objects
  return aliasDict

def get_genes_sizes(alias_dict, genomic_regions_file_name):

  # Initialization
  regionFile = open(genomic_regions_file_name, "rU")
  geneSizesDict = dict()

  # Iterating on region file
  for line in regionFile:

    # Fetching data
    ll = line.strip().split("\t")
    nn = ll[3].split(":")
    chrom = ll[0]; start = int(ll[1]); end = int(ll[2]); strand = ll[5]; region = nn[0]; name = nn[1]; activity = nn[2]
    if(region != "GENE" or activity == "INACTIVE"): continue

    # Gene name
    try: gene = alias_dict[name.upper()]
    except Exception: gene = name.upper()

    # Appending gene to allList
    geneSizesDict[gene] = end - start

  # Termination
  regionFile.close()

  # Returning objects
  return geneSizesDict

def get_promoters(promoter_ext, alias_dict, genomic_regions_file_name, stag_regions_file_name):

  # Initialization
  regionFile = open(genomic_regions_file_name, "rU")
  stagRegionsFile = Samfile(stag_regions_file_name, "rb")
  allList = []
  stagDict = dict()

  # Iterating on region file
  for line in regionFile:

    # Fetching data
    ll = line.strip().split("\t")
    nn = ll[3].split(":")
    chrom = ll[0]; start = ll[1]; end = ll[2]; strand = ll[5]; region = nn[0]; name = nn[1]; activity = nn[2]
    if(region != "PROMOTER" or activity == "INACTIVE"): continue

    # Gene name
    try: gene = alias_dict[name.upper()]
    except Exception: gene = name.upper()
    promoterWriteVec = [chrom, start, end, gene, "0", strand]

    # Appending gene to allList
    allList.append(promoterWriteVec)

    # Check whether promoter intersect both stag regions
    promoterRegion = [chrom, max(int(start)-promoter_ext, 0), int(end)+promoter_ext]
    check = check_bam_at_least_one_read(stagRegionsFile, promoterRegion)
    if(not check): continue
    stagDict[gene] = promoterWriteVec

  # Termination
  regionFile.close()
  stagRegionsFile.close()

  # Returning objects
  return allList, stagDict

def get_expression_dictionary(alias_dict, expression_file_name):
  
  # Fetching expression dictionary
  expression_file = open(expression_file_name, "rU")
  expression_file.readline()
  exprDict = dict()
  for line in expression_file:
    ll = line.strip().split("\t")
    try: gene = alias_dict[ll[1].upper()]
    except Exception: gene = ll[1].upper()
    exprDict[gene] = float(ll[2])
  expression_file.close()

  # Returning objects
  return exprDict

def get_stag_expression(stag_genes_dict, expression_dict):

  # Calculating STAG fold change expression dictionary
  fcDict = dict()
  for k in stag_genes_dict.keys():

    # Fetching gene name
    gene = stag_genes_dict[k][3]
    try: fc = expression_dict[gene]
    except Exception: continue
    fcDict[gene] = fc

  # Returning objects
  return fcDict

def get_control_expression(random_repeats, number_of_measures, all_promoters_list, expression_dict):

  # Fetching matrix with <random_repeats> x <number_of_measures> random expression fold-changes
  randomMatrix = []
  for i in range(0, random_repeats):
    randomVec = []
    while len(randomVec) < number_of_measures:
      r = randint(0, len(all_promoters_list) - 1)
      gene = all_promoters_list[r][3]
      
      # Fetching expression
      try: fc = expression_dict[gene]
      except Exception: continue
      randomVec.append(fc)

    # Appending vector
    randomMatrix.append(sorted(randomVec))

  # Summing the sorted vectors of the matrix
  resVec = [0.0] * len(randomMatrix[0])
  for col in range(0, len(resVec)):
    summ = 0.0
    for row in range(0, len(randomMatrix)): summ += randomMatrix[row][col]
    resVec[col] = summ / len(randomMatrix)

  # Returning objects
  return resVec

def write_table(stag_exp_dict, control_exp_list, output_table_file_name):

  # Writing table for R
  stagKeyList = sorted(stag_exp_dict.keys())
  header = ["GENE", "STAG_EXP", "CONT_EXP"]
  outputTableFile = open(output_table_file_name, "w")
  outputTableFile.write("\t".join(header)+"\n")
  for i in range(0, len(stagKeyList)):
    gene = stagKeyList[i]
    stagExp = stag_exp_dict[gene]
    contrExp = control_exp_list[i]
    outputTableFile.write("\t".join([str(e) for e in [gene, stagExp, contrExp]])+"\n")
  outputTableFile.close()

###################################################################################################
# Execution
###################################################################################################

# Fetch alias dictionary
aliasDict = get_alias_dictionary(aliasFileName)

# Fetch all genes
allPromotersList, stagGenesDict = get_promoters(promoterExt, aliasDict, genomicRegionsFileName, stagRegionsFileName)

# Fetch expression dictionary
expDict = get_expression_dictionary(aliasDict, expressionFileName)

# Fetch STAG expression
stagExpDict = get_stag_expression(stagGenesDict, expDict)

# Fetch CONTROL expression
controlExpList = get_control_expression(randomRepeats, len(stagExpDict), allPromotersList, expDict)

# Create table for R
write_table(stagExpDict, controlExpList, outputTableFileName)


