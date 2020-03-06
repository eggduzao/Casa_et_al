
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
halfBin = int(sys.argv[1])
resolution = int(sys.argv[2])
regionFileName = sys.argv[3]
matrix1FileName = sys.argv[4]
matrix2FileName = sys.argv[5]
outputTableFileName = sys.argv[6]

# Initialization
outputLoc = "/".join(outputTableFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def round_down(num, divisor):
  return num - (num%divisor)

def round_up(num, divisor):
  return num + (divisor - (num%divisor))

def read_region_list(half_bin, resolution, region_file_name):

  # Initialization
  region_list = []
  region_file = open(region_file_name, "rU")

  # Iterating in the input file
  for line in region_file:

    # Original region
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])
    mid = (p1+p2)/2
    
    # Upstream region
    r1 = round_down(mid, resolution) - (half_bin * resolution)

    # Downstream region
    r2 = round_up(mid, resolution) + (half_bin * resolution)

    # Appending matrix region to list
    region = [chrom, r1, r2]
    region_list.append(region)

  # Returning objects
  region_file.close()
  return region_list

def read_matrix_dictionary(matrix_file_name):

  # Initialization
  matrix_dict = dict()
  matrix_file = open(matrix_file_name,"rU")
  bad_list = ["NaN", "Inf", "-Inf", "NA", "0", "0.0"]

  # Populating matrix
  for line in matrix_file:
    ll = line.strip().split("\t")
    if(ll[3] in bad_list): continue
    p1 = int(ll[1])
    p2 = int(ll[2])
    minP = str(min(p1, p2))
    maxP = str(max(p1, p2))
    matrix_dict[":".join([ll[0], minP, maxP])] = float(ll[3])
    matrix_dict[":".join([ll[0], maxP, minP])] = float(ll[3])

  # Returning objects
  matrix_file.close()
  return matrix_dict

def update_matrix(half_bin, resolution, matrix, region, matrix_dict):

  # Initialization
  chrom = region[0]

  # Updating contact matrix with novel signal
  counterI = 0
  for i in range(region[1], region[2], resolution):
    counterJ = 0
    for j in range(region[1], region[2], resolution):
      try:
        value = matrix_dict[":".join([chrom, str(i), str(j)])]
        matrix[counterI][counterJ] = matrix[counterI][counterJ] + value
      except Exception: pass
      counterJ += 1
    counterI += 1

  # Returning status
  return 0

def update_matrix_fc(half_bin, resolution, matrix, region, matrix_1_dict, matrix_2_dict):

  # Initialization
  chrom = region[0]

  # Updating contact matrix with novel signal
  counterI = 0
  for i in range(region[1], region[2], resolution):
    counterJ = 0
    for j in range(region[1], region[2], resolution):
      flag_update_1 = True; flag_update_2 = True; value1 = 0; value2 = 0
      try: value1 = matrix_1_dict[":".join([chrom, str(i), str(j)])]
      except Exception: flag_update_1 = False
      try: value2 = matrix_2_dict[":".join([chrom, str(i), str(j)])]
      except Exception: flag_update_2 = False
      if(flag_update_1 or flag_update_2): matrix[counterI][counterJ] = matrix[counterI][counterJ] + (value1 - value2)
      counterJ += 1
    counterI += 1

  # Returning status
  return 0

def print_matrix(matrix, total_counter, output_table_file_name):

  # Printing final table containing the matrix
  output_table_file = open(output_table_file_name,"w")
  for vec in matrix:
    to_print = []
    for v in vec: to_print.append(str(v / total_counter))
    output_table_file.write("\t".join(to_print)+"\n")
  output_table_file.close()

  # Returning status
  return 0

def create_table(half_bin, resolution, region_file_name, matrix_1_file_name, matrix_2_file_name, output_table_file_name):

  # Reading regions
  region_list = read_region_list(half_bin, resolution, region_file_name) 

  # Reading matrix 1
  matrix_1_dict = read_matrix_dictionary(matrix_1_file_name)

  # Reading matrix 2
  fold_change_flag = False
  if(os.path.isfile(matrix_2_file_name)):
    fold_change_flag = True
    matrix_2_dict = read_matrix_dictionary(matrix_2_file_name)

  # Creating main matrix
  matrix = [[0.0] * ((half_bin*2)+1) for e in range(0,((half_bin*2)+1))]

  # Loop on regions to create table
  total_counter = 0.0
  for region in region_list:

    # Updating matrix
    if(fold_change_flag): update_matrix_fc(half_bin, resolution, matrix, region, matrix_1_dict, matrix_2_dict)
    else: update_matrix(half_bin, resolution, matrix, region, matrix_1_dict)

    # Update average counter
    total_counter += 1.0

  # Printing table
  print_matrix(matrix, total_counter, output_table_file_name) 

###################################################################################################
# Execution
###################################################################################################

# Create Table
create_table(halfBin, resolution, regionFileName, matrix1FileName, matrix2FileName, outputTableFileName)


