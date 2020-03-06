
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import math
sys.path = ["/home/egusmao/.local/lib/python2.7/site-packages"] + sys.path
import pyBigWig
from pysam import Samfile
from numpy import mean, array_split

# Input
downstream_extension = int(sys.argv[1])
upstream_extension = int(sys.argv[2])
numberOfBins = int(sys.argv[3])
numberOfCounts = int(sys.argv[4])
aliasFileName = sys.argv[5]
geneListFileName = sys.argv[6]
regionsBedFileName = sys.argv[7]
signalFileType = sys.argv[8]
signalFileName = sys.argv[9]
outputGenesFileName = sys.argv[10]

###################################################################################################
# Functions
###################################################################################################

def create_alias_dictionary(alias_file_name):

  # Creating alias dictionary
  aliasDict = dict()
  aliasFile = open(alias_file_name,"rU")
  for line in aliasFile:
    ll = line.strip().split("\t")
    value = ll[1]
    geneList = [ll[0],ll[1]]+ll[2].split("&")
    for g in geneList: aliasDict[g.upper()] = value.upper()
  aliasFile.close()

  # Returning objects
  return aliasDict

def create_gene_list_dictionary(alias_dict, gene_list_file_name):

  # Structures
  gene_list_dict = dict() # GENE SYMBOL -> TRUE
  
  # Creating gene list dictionary
  if(gene_list_file_name):
    gene_list_file = open(gene_list_file_name, "rU")
    for line in gene_list_file:
      gene = line.strip().upper()
      try: gene = alias_dict[gene]
      except Exception: pass
      gene_list_dict[gene] = True
    gene_list_file.close()
  
    # Returning objects
    return gene_list_dict

  else:

    # Returning objects
    return None

def genes_dictionary(alias_dict, gene_list_dict, regions_file_name):

  # Structures
  # geneDict = GENE_NAME -> [CHROMOSOME, START, END, GENE, <ACTIVE (1) | INACTIVE (0)>, STRAND]
  geneDict = dict()

  # Iterating over the regions file populating the dictionaries
  regions_file = open(regions_file_name, "rU")
  for line in regions_file:
    ll = line.strip().split("\t")
    nn = ll[3].split(":")
    chrom = ll[0]; start = ll[1]; end = ll[2]; strand = ll[5]
    region = nn[0]; name = nn[1].upper(); activity = nn[2]
    if(region != "GENE"): continue
    try: gene = alias_dict[name]
    except Exception: gene = name
    if(gene_list_dict):
      try: is_the_gene = gene_list_dict[gene]
      except Exception: continue
    if(activity == "ACTIVE"): activity = "1"
    else: activity = "0"
    geneDict[gene] = [chrom, start, end, gene, activity, strand]
  regions_file.close()

  # Returning objects
  return geneDict

def fetch_counts(downstream_extension, upstream_extension, value_to_add, input_file, region, filetype = "bam"):

  # Creating vector
  vector = []

  if(filetype == "bam"):

    # Creating vector
    total_length = (region[2]-region[1])
    incr = 2 * (downstream_extension + upstream_extension)
    vector = [0.0] * (total_length+(2*incr))
    region1i = region[1] - incr
    region2i = region[2] + incr
    if(region1i < 0): return None

    # Fetching bam signal
    for read in input_file.fetch(region[0], region1i, region2i):
      if(read.is_reverse):
        read_start = max(read.reference_end - downstream_extension, region1i)
        read_end = min(read.reference_end + upstream_extension, region2i)
      else:
        read_start = max(read.reference_start - upstream_extension, region1i)
        read_end = min(read.reference_start + downstream_extension, region2i)
      for i in range(read_start, read_end): vector[i-region1i] += value_to_add

    # Trimming vector
    vector = vector[(incr):(len(vector)-2*incr)]
    for i in range(0, len(vector)):
      if(not vector or math.isnan(vector[i]) or math.isinf(vector[i])): vector[i] = 0.0

  elif(filetype == "bw"):

    # Creating vector
    total_length = (region[2]-region[1])
    vector = [0.0] * total_length

    # Fetching bw signal
    try: valuesVec = input_file.values(region[0], region[1], region[2])
    except Exception: return vector
    for i in range(0, len(valuesVec)):
      if(not valuesVec[i] or math.isnan(valuesVec[i]) or math.isinf(valuesVec[i])): vector[i] = 0.0
      else: vector[i] = valuesVec[i]

  # Returning objects
  return vector

def fetch_binarized_vector(number_of_bins, signal_vector):

  # Returning objetcs
  try:
    return [mean(e) if sum(e) > 0 else 0.0 for e in array_split(signal_vector, number_of_bins)]
  except Exception:
    return None

def runing_average(vector, smoothing_number):

  # Passing running means
  newVector = []
  for i in range(0, len(vector)): newVector.append(mean(vector[max(0,i-smoothing_number):min(len(vector),i+smoothing_number)]))

  # Returning objetcs
  return newVector

def write_meta_signals(downstream_extension, upstream_extension, number_of_bins, number_of_counts, region_dict, signal_file_type, signal_file, output_file_name):

  # Calculating rpm and valueToAdd
  rpm = float(number_of_counts) / 1000000.
  valueToAdd = 1.0 / rpm

  # Opening files
  output_file = open(output_file_name, "w")

  # Iterating on dictionary
  for rkey in region_dict.keys():

    # Coordinates and status
    r = region_dict[rkey]
    chrom = r[0]; start = int(r[1]); end = int(r[2]); gene = r[3]; activity = r[4]; strand = r[5]

    # Binning based on gene
    step = 5000
    tss_2 = max(start - (2*step), 0); tss_1 = max(start - step, 0); tss = start; tes = end; tes_1 = end + step; tes_2 = end + (2*step)

    # Fetch signal on each bin
    signal1 = fetch_binarized_vector(number_of_bins, fetch_counts(downstream_extension, upstream_extension, valueToAdd, signal_file, [chrom, tss_2, tss_1], filetype = signal_file_type))
    signal2 = fetch_binarized_vector(number_of_bins, fetch_counts(downstream_extension, upstream_extension, valueToAdd, signal_file, [chrom, tss_1, tss], filetype = signal_file_type))
    signal3 = fetch_binarized_vector(number_of_bins, fetch_counts(downstream_extension, upstream_extension, valueToAdd, signal_file, [chrom, tss, tes], filetype = signal_file_type))
    signal4 = fetch_binarized_vector(number_of_bins, fetch_counts(downstream_extension, upstream_extension, valueToAdd, signal_file, [chrom, tes, tes_1], filetype = signal_file_type))
    signal5 = fetch_binarized_vector(number_of_bins, fetch_counts(downstream_extension, upstream_extension, valueToAdd, signal_file, [chrom, tes_1, tes_2], filetype = signal_file_type))
    if(signal1 and signal2 and signal3 and signal4 and signal5):
      signal = signal1 + signal2 + signal3 + signal4 + signal5
    else: continue
    if(strand == "-"): signal = signal[::-1]
    final_signal = runing_average(signal, step/numberOfBins)

    # Writing output
    start_values = [chrom, start, end, gene, activity, strand]
    output_file.write("\t".join([str(e) for e in start_values + signal])+"\n")

  # Closing file
  output_file.close()

def create_table(downstream_extension, upstream_extension, number_of_bins, number_of_counts, alias_file_name, gene_list_file_name, regions_file_name, signal_file_type, signal_file_name, output_gene_file_name):

  # Fetch alias dictionary
  aliasDict = create_alias_dictionary(alias_file_name)

  # Fetch gene list dict
  if(gene_list_file_name == "."): gene_list_file_name = None
  geneListDict = create_gene_list_dictionary(aliasDict, gene_list_file_name)

  # Fetch genes and enhancers
  geneDict = genes_dictionary(aliasDict, geneListDict, regions_file_name)

  # Opening signal file
  if(signal_file_type == "bam"): 
    signal_file = Samfile(signal_file_name, "rb")
  elif(signal_file_type == "bw"):
    signal_file = pyBigWig.open(signal_file_name)

  # Writing meta signals
  write_meta_signals(downstream_extension, upstream_extension, number_of_bins, number_of_counts, geneDict, signal_file_type, signal_file, output_gene_file_name)
  signal_file.close()

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
create_table(downstream_extension, upstream_extension, numberOfBins, numberOfCounts, aliasFileName, geneListFileName, regionsBedFileName, signalFileType, signalFileName, outputGenesFileName)


