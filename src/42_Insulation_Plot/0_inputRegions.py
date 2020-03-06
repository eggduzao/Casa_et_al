
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
region = sys.argv[1]
aliasDictFileName = sys.argv[2]
regionsFileName = sys.argv[3]
genesFileName = sys.argv[4]
tempLocation = sys.argv[5]
outputFileName = sys.argv[6]

# Initialization
promExt = 2000
outputLoc = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLoc
os.system(command)
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def sort(input_file, output_file):
  command = "sort -k1,1 -k2,2n "+input_file+" > "+output_file
  os.system(command)

def create_alias_dictionary(alias_file_name):

  # Create structures
  aliasDict = dict()
  aliasFile = open(alias_file_name,"rU")

  # Iterate over file to populate dictionary
  for line in aliasFile:
    ll = line.strip().split("\t")
    value = ll[1]
    geneList = [ll[0],ll[1]]+ll[2].split("&")
    for g in geneList: aliasDict[g.upper()] = value.upper()

  # Returning objects
  aliasFile.close()
  return aliasDict

def create_gene_dictionary(alias_dict, genes_file_name):

  # Create structures
  genes_dict = dict() # Gene -> True
  genes_file = open(genes_file_name,"rU")

  # Iterate over file to populate dictionary
  for line in genes_file:
    geneTemp = line.strip().upper()
    try: gene = alias_dict[geneTemp]
    except Exception: gene = geneTemp
    genes_dict[gene] = True

  # Returning objects
  genes_file.close()
  return genes_dict

def create_promoter_list(region, alias_dict, genes_dict, regions_file_name):

  # Create structures
  region_list = []
  regions_file = open(regions_file_name, "rU")

  # Iterate over file to get all TSSs
  for line in regions_file:
    ll = line.strip().split("\t")
    nn = ll[3].split(":")
    chrom = ll[0]; start = ll[1]; end = ll[2]; strand = ll[5]
    region = nn[0]; name = nn[1].upper(); activity = nn[2]
    if(region != "GENE"): continue
    try: gene = alias_dict[name]
    except Exception: gene = name
    try: keep = genes_dict[gene]
    except Exception: continue
    if(region == "TSS"):
      if(strand == "+"): region_list.append([chrom, str(max(int(start) - promExt, 0)), start, gene, "1", strand])
      else: region_list.append([chrom, end, str(int(end) + promExt), gene, "1", strand])
    else:
      if(strand == "+"): region_list.append([chrom, end, str(int(end) + promExt), gene, "2", strand])
      else: region_list.append([chrom, str(max(int(start) - promExt, 0)), start, gene, "2", strand])

  # Returning objects
  regions_file.close()
  return region_list

###################################################################################################
# Execution
###################################################################################################

# Read alias dictionary
alias_dict = create_alias_dictionary(aliasDictFileName)

# Read genes dictionary
genes_dict = create_gene_dictionary(alias_dict, genesFileName)

# Read promoter list
promoter_list = create_promoter_list(region, alias_dict, genes_dict, regionsFileName)

# Opening output temporary files
outputFileTempName = tempLocation + "outputFileTempName.bed"
outputFileTemp = open(outputFileTempName, "w")

# Iterating on promoters and writing regions
for r in promoter_list:

  # Initialization
  chrom = r[0]; start = r[1]; end = r[2]; gene = r[3]; strand = r[5]

  # Writing region
  outputFileTemp.write("\t".join([chrom, start, end, gene, "1", strand])+"\n")

# Closing files
outputFileTemp.close()

# Sorting files
sort(outputFileTempName, outputFileName)


