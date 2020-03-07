
# Import
import os
import sys

###################################################################################################
# INPUT
###################################################################################################

# Input
genomeId = sys.argv[1]
chromosome = sys.argv[2]
resolution = sys.argv[3]
juicerCommand = sys.argv[4]
inputSparseFileName = sys.argv[5]
temporaryLocation = sys.argv[6]
outputHicFileName = sys.argv[7]

# Initialization
command = "mkdir -p "+temporaryLocation
os.system(command)
outputLocation = "/".join(outputHicFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# FUNCTIONS
###################################################################################################

def read_sparse_matrix(input_sparse_file_name):
  
  # Reading matrix
  sparse_matrix_dict = dict()
  input_sparse_file = open(input_sparse_file_name, "rU")
  for line in input_sparse_file:
    ll = line.strip().split("\t")
    key = ":".join(ll[:3])
    sparse_matrix_dict[key] = ll[3]
  input_sparse_file.close()

  # Returning objects
  return sparse_matrix_dict

def sparse_matrix_to_pre(sparse_matrix_dictionary, pre_file_name):

  # Writing matrix to pre_file_name
  pre_file = open(pre_file_name, "w")
  for key in sparse_matrix_dictionary.keys():
    kk = key.split(":")
    str1 = "0"; str2 = "1"; frag1 = "0"; frag2 = "1"
    chr1 = kk[0]; chr2 = kk[0]; pos1 = kk[1]; pos2 = kk[2]
    score = sparse_matrix_dictionary[key]
    vector = [str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, score]
    pre_file.write(" ".join(vector)+"\n")

  # Termination
  pre_file.close()

def sort_pre_file(pre_file_name, temporary_location, pre_file_name_sorted):

  # 1 - for the first read end chromosome to be less than the second read end chromosome
  temp_sorted_file_name = temporary_location + "temp_sorted_file_name.pre"
  command = "awk '{if ($3 > $7){ print $5, $6, $7, $8, $1, $2, $3, $4, $9}else {print}}' " + pre_file_name + " > " + temp_sorted_file_name
  os.system(command)

  # 2 - for the reads to be sorted by chromosome block. That is, all chr3R-chr3R reads together in one place. This is so we don’t have to read the file multiple times.
  command = "sort -k2,2d -k6,6d " + temp_sorted_file_name + " > " + pre_file_name_sorted
  os.system(command)

def pre_to_hic(genome_id, chromosome, resolution, juicer_command, pre_file_name, temporary_location, output_hic_file_name):

  # Apply pre
  command = " ".join([juicer_command, "pre", "-d", "-n", "-r", resolution, "-c", chromosome, "-t", temporary_location, pre_file_name, output_hic_file_name, genome_id])
  os.system(command)

def merge_single_cell_files(genome_id, chromosome, resolution, juicer_command, input_sparse_file_name, temporary_location, output_hic_file_name):

  # Get matrix dictionary
  sparse_matrix_dict = read_sparse_matrix(input_sparse_file_name)

  # Convert to unsorted pre file
  unsorted_pre_file_name = temporary_location + "unsorted_pre_file_name.pre"
  sparse_matrix_to_pre(sparse_matrix_dict, unsorted_pre_file_name)

  # Sort pre file name
  sorted_pre_file_name = temporary_location + "sorted_pre_file_name.pre"
  sort_pre_file(unsorted_pre_file_name, temporary_location, sorted_pre_file_name)

  # Convert pre file to hic
  pre_to_hic(genome_id, chromosome, resolution, juicer_command, sorted_pre_file_name, temporary_location, output_hic_file_name)

###################################################################################################
# EXECUTION
###################################################################################################

merge_single_cell_files(genomeId, chromosome, resolution, juicerCommand, inputSparseFileName, temporaryLocation, outputHicFileName)


