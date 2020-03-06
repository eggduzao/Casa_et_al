
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

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

def only(input_file1, input_file2, output_file):
  command = "intersectBed -wa -v -a "+input_file1+" -b "+input_file2+" > "+output_file
  os.system(command)

###################################################################################################
# Creating file list
###################################################################################################

# Parameters
ll = "/projects/ag-papan/eduardo/Wendt_Stag/Results/9_Process_Losada_Data/4_Peaks/"
aliasFileName = "/home/egusmao/rgtdata/hg19/alias_human_booster.txt"
tempLoc = "./TEMP/"
regionsBedFileName = "/projects/ag-papan/eduardo/Wendt_Stag/Results/0_Definitive_Gene_Annotation/regions.bed"
losadaStag1FileName = ll + "SA1_MCF10A_CHIPSEQ_peaks_filtered.narrowPeak"
losadaStag2FileName = ll + "SA2_MCF10A_CHIPSEQ_peaks_filtered.narrowPeak"
il = "/projects/ag-papan/eduardo/Wendt_Stag/Results/10_Losada_Meta_Plots_Genes/0_stag_bed_files_losada/"

# Parameters
command = "mkdir -p "+tempLoc
os.system(command)

# Cut
losadaStag1FileNameCut = tempLoc + "losadaStag1FileNameCut.bed"
losadaStag2FileNameCut = tempLoc + "losadaStag2FileNameCut.bed"
cut(losadaStag1FileName, losadaStag1FileNameCut)
cut(losadaStag2FileName, losadaStag2FileNameCut)

# Full Peaks 
losadaStag1FileNameSort = il + "stag1_all.bed"
losadaStag2FileNameSort = il + "stag2_all.bed"
cut(losadaStag1FileNameCut, losadaStag1FileNameSort)
cut(losadaStag2FileNameCut, losadaStag2FileNameSort)

# Shared
sharedFileName = il + "shared.bed"
intersect(losadaStag1FileNameSort, losadaStag2FileNameSort, sharedFileName)

# Only
stag1FileName = il + "stag1_only.bed"
stag2FileName = il + "stag2_only.bed"
only(losadaStag1FileNameSort, losadaStag2FileNameSort, stag1FileName)
only(losadaStag2FileNameSort, losadaStag1FileNameSort, stag2FileName)

###################################################################################################
# Execution
###################################################################################################

# All files
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Data/losada/expression/gene_lists/"
stagListFileName = [losadaStag1FileNameSort, losadaStag2FileNameSort, sharedFileName, stag1FileName, stag2FileName]
labelList = ["stag1_all", "stag2_all", "shared", "stag1_only", "stag2_only"]

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
  stagFileName = stagListFileName[i]
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


