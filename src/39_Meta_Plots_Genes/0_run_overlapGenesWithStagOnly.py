
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Stag Regions List
il = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/stag_bed_files/"
ol = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/expression/gene_lists/"
aliasFileName = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
tempLoc = "./TEMP/"
stagListFileName = ["nonpredominant", "shared", "STAG1_full_peaks", "STAG1_only", "STAG1_predominant", "STAG2_full_peaks", "STAG2_only", "STAG2_predominant"]
labelList = ["nonpredominant", "shared", "stag1_all", "stag1_only", "stag1_predominant", "stag2_all", "stag2_only", "stag2_predominant"]
regionsBedFileName = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/0_Definitive_Gene_Annotation/regions.bed"

# Parameters
command = "mkdir -p "+tempLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def create_alias_dictionary(alias_file_name):
  aliasDict = dict()
  aliasFile = open(alias_file_name,"rU")
  for line in aliasFile:
    ll = line.strip().split("\t")
    value = ll[1]
    geneList = [ll[0],ll[1]]+ll[2].split("&")
    for g in geneList: aliasDict[g.upper()] = value.upper()
  aliasFile.close()
  return aliasDict

def create_genes_only_file(alias_dict, regions_file_name, output_file_name):
  regions_file = open(regions_file_name, "rU")
  output_file = open(output_file_name, "w")
  for line in regions_file:
    ll = line.strip().split("\t")
    nn = ll[3].split(":")
    chrom = ll[0]; start = ll[1]; end = ll[2]; strand = ll[5]
    region = nn[0]; name = nn[1].upper(); activity = nn[2]
    if(region != "GENE"): continue
    try: gene = alias_dict[name]
    except Exception: gene = name
    output_file.write("\t".join([chrom, start, end, gene, "1", strand])+"\n")
  regions_file.close()
  output_file.close()

def cut(input_file, output_file):
  command = "cut -f 1,2,3 "+input_file+" > "+output_file
  os.system(command)

def sort(input_file, output_file):
  command = "sort -k1,1 -k2,2n "+input_file+" > "+output_file
  os.system(command)

def intersect(input_file1, input_file2, output_file):
  command = "intersectBed -wa -u -a "+input_file1+" -b "+input_file2+" > "+output_file
  os.system(command)

###################################################################################################
# Execution
###################################################################################################

# Creating alias dictionary
aliasDict = create_alias_dictionary(aliasFileName)

# Create genes file
genesFileName = tempLoc + "genesFileName.bed"
create_genes_only_file(aliasDict, regionsBedFileName, genesFileName)
genesSortedFileName = tempLoc + "genesSortedFileName.bed"
sort(genesFileName, genesSortedFileName)

# Stag Regions Loop
for i in range(0, len(stagListFileName)):

  # Input
  stagFileName = il + stagListFileName[i] + ".bed"
  outputFileName = ol + labelList[i] + ".txt"

  # Correcting stag file
  cutStagFileName = tempLoc + "cutStagFileName.bed"
  cut(stagFileName, cutStagFileName)
  sortStagFileName = tempLoc + "sortStagFileName.bed"
  sort(cutStagFileName, sortStagFileName)

  # Intersecting stag and genes
  intersectFileName = tempLoc + "intersectFileName.bed"
  intersect(genesSortedFileName, sortStagFileName, intersectFileName)

  # Writing genes to list
  intersectFile = open(intersectFileName, "rU")
  outputFile = open(outputFileName, "w")
  for line in intersectFile: outputFile.write(line.strip().split("\t")[3]+"\n")
  intersectFile.close()
  outputFile.close()

# Removing tempLoc
command = "rm -rf "+tempLoc
os.system(command)


