
# Import
import os
import sys
from math import isnan, isinf

# Input
genome_id = sys.argv[1]
resolution = sys.argv[2]
chrom_sizes_file_name = sys.argv[3]
input_hic_file_name_list = sys.argv[4].split(",")
temporary_location = sys.argv[5]
output_hic_file_name = sys.argv[6]
output_sparse_file_name = sys.argv[7]

# Parameters
juicer_command = "juicertools"
kind_of_matrix = "observed"
kind_of_normalization = "NONE"
unit_of_resolution = "BP"

# Initialization
command = "mkdir -p "+temporary_location
os.system(command)
outputLocation = "/".join(output_hic_file_name.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def sort_pre_file(pre_file_name, temporary_location, pre_file_name_sorted):

  # 1 - for the first read end chromosome to be less than the second read end chromosome
  temp_sorted_file_name = temporary_location + "temp_sorted_file_name.pre"
  command = "awk '{if ($3 > $7){ print $5, $6, $7, $8, $1, $2, $3, $4, $9}else {print}}' " + pre_file_name + " > " + temp_sorted_file_name
  os.system(command)

  # 2 - for the reads to be sorted by chromosome block. That is, all chr3R-chr3R reads together in one place. This is so we don’t have to read the file multiple times.
  command = "sort -k2,2d -k6,6d " + temp_sorted_file_name + " > " + pre_file_name_sorted
  os.system(command)

def get_chromosome_sizes_dictionary(chrom_sizes_file_name):

  # Chromosome sizes dict
  chrom_size_dict = dict()
  chrom_sizes_file = open(chrom_sizes_file_name,"rU")
  for line in chrom_sizes_file:
    ll = line.strip().split("\t")
    chrom_size_dict[ll[0]] = int(ll[1])
  chrom_sizes_file.close()
  chrom_list = sorted(chrom_size_dict.keys())

  # Return objects
  return chrom_size_dict, chrom_list

def hic_to_sparse_matrix(juicer_command, kind_of_matrix, kind_of_normalization, unit_of_resolution, resolution, chromosome_list, chrom_size_dict, input_hic_file_name, temporary_location, sparse_matrix_dictionary):

  # Creating individual sparse chromosome matrices
  for chrom in chromosome_list:

    # Initialization
    chrWoChr = chrom.split("chr")[-1]
    region = ":".join([chrWoChr,"1",str(chrom_size_dict[chrom])])
    output_chrom_file_name = temporary_location + chrom + ".txt"

    # Creating sparse matrix
    command = " ".join([juicer_command, "dump", kind_of_matrix, kind_of_normalization, input_hic_file_name, region, region, unit_of_resolution, resolution, output_chrom_file_name])
    os.system(command)

    # Writing entries to sparse_matrix_dictionary
    output_chrom_file = open(output_chrom_file_name, "rU")
    for line in output_chrom_file:
      ll = line.strip().split("\t")
      try: value = float(ll[2])
      except Exception: continue
      if(isnan(value) or isinf(value)): continue
      key = ":".join([chrom, ll[0], ll[1]])
      try: sparse_matrix_dictionary[key] += value
      except Exception: sparse_matrix_dictionary[key] = value

def sparse_matrix_to_pre(sparse_matrix_dictionary, pre_file_name, output_sparse_file_name):

  # Opening sparse matrix file
  output_sparse_file = open(output_sparse_file_name, "w")

  # Writing matrix to pre_file_name
  pre_file = open(pre_file_name, "w")
  for key in sparse_matrix_dictionary.keys():
    kk = key.split(":")
    str1 = "0"; str2 = "1"
    frag1 = "0"; frag2 = "1"
    chr1 = kk[0]; chr2 = kk[0]
    if(chr1 != chr2 or chr1 == "M" or chr1 == "Y" or chr2 == "M" or chr2 == "Y"): continue
    pos1 = kk[1]; pos2 = kk[2]
    score = str(sparse_matrix_dictionary[key])
    vector = [str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, score]
    pre_file.write(" ".join(vector)+"\n")
    vectorsp = [chr1, pos1, pos2, score]
    output_sparse_file.write("\t".join(vectorsp)+"\n")

  # Termination
  pre_file.close()
  output_sparse_file.close()

def pre_to_hic(genome_id, resolution, juicer_command, pre_file_name, temporary_location, output_hic_file_name):

  # Apply pre
  command = " ".join([juicer_command, "pre", "-d", "-n", "-r", resolution, "-t", temporary_location, pre_file_name, output_hic_file_name, genome_id])
  os.system(command)

def merge_single_cell_files(genome_id, juicer_command, kind_of_matrix, kind_of_normalization, unit_of_resolution, resolution, chrom_sizes_file_name, input_hic_file_name_list, temporary_location, output_hic_file_name, output_sparse_file_name):

  # Fetching chromosome sizes
  chrom_size_dict, chromosome_list = get_chromosome_sizes_dictionary(chrom_sizes_file_name)

  # Creating sparse matrix
  sparse_matrix_dictionary = dict()

  # Iterating on the input_hic_file_name_list
  for i in range(0, len(input_hic_file_name_list)):

    # Parameters
    input_hic_file_name = input_hic_file_name_list[i]

    # Adding values to the sparse matrix
    hic_to_sparse_matrix(juicer_command, kind_of_matrix, kind_of_normalization, unit_of_resolution, resolution, chromosome_list, chrom_size_dict, input_hic_file_name, temporary_location, sparse_matrix_dictionary)

  # Creating pre file
  pre_file_name = temporary_location + "pre_file_name.pre"
  sparse_matrix_to_pre(sparse_matrix_dictionary, pre_file_name, output_sparse_file_name)

  # Sorting pre file
  pre_file_name_sorted = temporary_location + "pre_file_name_sorted.pre"
  sort_pre_file(pre_file_name, temporary_location, pre_file_name_sorted)

  # Creating final hic file
  pre_to_hic(genome_id, resolution, juicer_command, pre_file_name_sorted, temporary_location, output_hic_file_name)

###################################################################################################
# EXECUTION
###################################################################################################

# Merge hic files
merge_single_cell_files(genome_id, juicer_command, kind_of_matrix, kind_of_normalization, unit_of_resolution, resolution, chrom_sizes_file_name, input_hic_file_name_list, temporary_location, output_hic_file_name, output_sparse_file_name)


