
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import gc

# Input
resolution = int(sys.argv[1])
percentileThreshold = int(sys.argv[2])
percentileFileName = sys.argv[3]
inputMatrixFileName = sys.argv[4]
tempLocation = sys.argv[5]
outBinPrefix = sys.argv[6]
outStdPrefix = sys.argv[7]

# Initialization
genomeId = "hg19"
command = "mkdir -p "+tempLocation
os.system(command)
outputLocation = "/".join(outBinPrefix.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def apply_pre(resolution, genome_id, pre_reads_file_name, temporary_location, output_hic_file_name):

  # Apply pre
  # -n Don't normalize the matrices
  # -d Only calculate intra chromosome
  # -r <comma-separated list of resolutions> Only calculate specific resolutions
  # -t <tmpDir> Set a temporary directory for writing
  command = "juicertools pre -n -d -r "+str(resolution)+" -t "+temporary_location+" "+" ".join([pre_reads_file_name,
            output_hic_file_name, genome_id])
  os.system(command)

def read_percentile_file(perc_file_name):

  perc_dict = dict()
  perc_file = open(perc_file_name, "rU")
  for line in perc_file:
    ll = line.strip().split("\t")
    perc_dict[int(ll[0])] = float(ll[1])
  perc_file.close()

  return perc_dict

def create_matrices(resolution, perc_threshold, perc_file_name, genome_id, temporary_location, input_matrix_file_name, out_bin_prefix, out_std_prefix):

  # File Names
  pre_bin_file_name = temporary_location + "pre_bin_file_name.txt"
  pre_std_file_name = temporary_location + "pre_std_file_name.txt"
  out_bin_text_file_name = out_bin_prefix + ".txt"
  out_std_text_file_name = out_std_prefix + ".txt"
  out_bin_hic_file_name = out_bin_prefix + ".hic"
  out_std_hic_file_name = out_std_prefix + ".hic"

  # Reading percentile file
  perc_dict = read_percentile_file(perc_file_name)

  # Opening all files
  input_matrix_file = open(input_matrix_file_name, "rU")
  pre_bin_file = open(pre_bin_file_name, "w")
  pre_std_file = open(pre_std_file_name, "w")
  out_bin_text_file = open(out_bin_text_file_name, "w")
  out_std_text_file = open(out_std_text_file_name, "w")

  # Verify minimum and maximum values of matrix outside diagonal to make standardization
  minV = 999999999.9
  maxV = -999999999.9

  # First iteration on matrix
  # 1 - Calculating minimum and maximum values outside diagonal
  # 2 - Writing binarized matrix
  # 3 - Creating binarized "pre" input
  for line in input_matrix_file:
   
    # Initialization
    ll = line.strip().split("\t")
    chrom = ll[0]; pos1 = ll[1]; pos2 = ll[2]; value = float(ll[3])

    # Evaluating minimum and maximum
    if(((int(pos1) < (int(pos2) - (5*resolution))) or (int(pos1) > (int(pos2) - (5*resolution)))) and (value > perc_dict[perc_threshold])):
      if(value < minV): minV = value
      if(value > maxV): maxV = value

    # Binarized matrix
    binValue = 0.0
    if(value > perc_dict[perc_threshold]): binValue = 1.0
    out_bin_text_file.write("\t".join([chrom, pos1, pos2, str(binValue)])+"\n")

    # Binarized "pre" input
    pre_bin_file.write(" ".join(["0", chrom, pos1, "0", "0", chrom, pos2, "1", str(binValue)])+"\n")

  # Closing files and cleaning environment
  out_bin_text_file.close()
  pre_bin_file.close()
  input_matrix_file.close()
  out_bin_text_file = None
  pre_bin_file = None
  input_matrix_file = None
  gc.collect()

  # Opening input file again
  input_matrix_file = open(input_matrix_file_name, "rU")

  # Second iteration on matrix
  # 1 - Writing standardized matrix
  # 2 - Creating standardized "pre" input
  for line in input_matrix_file:
   
    # Initialization
    ll = line.strip().split("\t")
    chrom = ll[0]; pos1 = ll[1]; pos2 = ll[2]; value = float(ll[3])  

    # Standardized matrix
    if(((int(pos1) < (int(pos2) - (5*resolution))) or (int(pos1) > (int(pos2) - (5*resolution)))) and (value > perc_dict[perc_threshold])):
      stdValue = (value - minV) / (maxV - minV)
      out_std_text_file.write("\t".join([chrom, pos1, pos2, str(stdValue)])+"\n")

      # Standardized "pre" input
      pre_std_file.write(" ".join(["0", chrom, pos1, "0", "0", chrom, pos2, "1", str(stdValue)])+"\n")

  # Closing files and cleaning environment
  out_std_text_file.close()
  pre_std_file.close()
  input_matrix_file.close()
  out_std_text_file = None
  pre_std_file = None
  input_matrix_file = None
  gc.collect()

  # Creating hic files
  apply_pre(resolution, genome_id, pre_bin_file_name, temporary_location, out_bin_hic_file_name)
  apply_pre(resolution, genome_id, pre_std_file_name, temporary_location, out_std_hic_file_name)

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
create_matrices(resolution, percentileThreshold, percentileFileName, genomeId, tempLocation, inputMatrixFileName, outBinPrefix, outStdPrefix)


