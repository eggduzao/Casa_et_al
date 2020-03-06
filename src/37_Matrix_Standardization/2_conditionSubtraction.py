
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import gc

# Input
resolution = sys.argv[1]
inputMatrix1FileName = sys.argv[2]
inputMatrix2FileName = sys.argv[3]
tempLocation = sys.argv[4]
outputMatrixPrefix = sys.argv[5]

# Initialization
genomeId = "hg19"
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def apply_pre(resolution, genome_id, pre_reads_file_name, temporary_location, output_hic_file_name):

  # Creating temporary sorted file
  out_name = output_hic_file_name.split("/")[-1].split(".")[0]
  temp_sorted_file_name = temporary_location + out_name + "_sorted.txt"
  command = "LC_ALL=C sort -k2,2 -k6,6 -k3,3g -k7,7g "+pre_reads_file_name+" > "+temp_sorted_file_name
  os.system(command)

  # Apply pre
  # -n Don't normalize the matrices
  # -d Only calculate intra chromosome
  # -r <comma-separated list of resolutions> Only calculate specific resolutions
  # -t <tmpDir> Set a temporary directory for writing
  command = "juicertools pre -n -d -r "+resolution+" -t "+temporary_location+" "+" ".join([temp_sorted_file_name,
            output_hic_file_name, genome_id])
  os.system(command)

def hic_signal_dictionary(input_matrix_file_name):
  
  # Execution
  resDict = dict() # CHROM:POS1:POS2 -> VALUE
  input_matrix_file = open(input_matrix_file_name, "rU")
  for line in input_matrix_file:
    ll = line.strip().split("\t")
    key = ":".join(ll[:3])
    resDict[key] = float(ll[3])
  input_matrix_file.close()
  return resDict

def subtract_matrices(resolution, genome_id, temporary_location, input_matrix_1_file_name, input_matrix_2_file_name, output_matrix_prefix):

  # File Names
  pre_file_name = temporary_location + "pre_file_name.txt"
  out_text_file_name = output_matrix_prefix + ".txt"
  out_hic_file_name = output_matrix_prefix + ".hic"

  # Opening all files
  pre_file = open(pre_file_name, "w")
  out_text_file = open(out_text_file_name, "w")

  # Fetching matrix signal dictionaries
  hicSignalDict1 = hic_signal_dictionary(input_matrix_1_file_name)
  hicSignalDict2 = hic_signal_dictionary(input_matrix_2_file_name)

  # Fetching list of all dictionaries' keys
  dictKeys = list(set(hicSignalDict1.keys())|set(hicSignalDict2.keys()))
  
  # Performing subtraction
  for key in dictKeys:
   
    # Initialization
    kk = key.split(":")
    chrom = kk[0]; pos1 = kk[1]; pos2 = kk[2]

    # Fetching values
    try: value1 = hicSignalDict1[key]
    except Exception: value1 = 0.0
    try: value2 = hicSignalDict2[key]
    except Exception: value2 = 0.0
    subValue = value1 - value2
    if(subValue == 0.0): continue

    # Writing text matrix
    out_text_file.write("\t".join([chrom, pos1, pos2, str(subValue)])+"\n")

    # Writing "pre" input matrix
    pre_file.write(" ".join(["0", chrom, pos1, "0", "0", chrom, pos2, "1", str(round(float(subValue),6))])+"\n")

  # Closing files and cleaning environment
  out_text_file.close()
  pre_file.close()
  out_text_file = None
  pre_file = None
  hicSignalDict1 = None
  hicSignalDict2 = None
  dictKeys = None
  gc.collect()

  # Creating hic files
  apply_pre(resolution, genome_id, pre_file_name, temporary_location, out_hic_file_name)

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
subtract_matrices(resolution, genomeId, tempLocation, inputMatrix1FileName, inputMatrix2FileName, outputMatrixPrefix)


