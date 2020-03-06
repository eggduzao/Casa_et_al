
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from math import floor, ceil
from collections import OrderedDict
sys.path = ["/home/egusmao/.local/lib/python2.7/site-packages"] + sys.path
from pysam import Samfile
from numpy import log2, mean
import gc

# Input
resolution = int(sys.argv[1])
number_of_bins = int(sys.argv[2]) # Number of <resolution> bins to "meta" the TADs
tad_1_file_name = sys.argv[3] # X-axis. Bed file.
tad_2_file_name = sys.argv[4] # Y-axis. Bed file.
hic_1_file_name = sys.argv[5]
hic_2_file_name = sys.argv[6] # If given, will calculate hic1 - hic2. If not, only hic1.
output_file_name = sys.argv[7]

###################################################################################################
# Functions
###################################################################################################

def fetch_tad_list(resolution, tad_file_name):

  # Read TAD list
  tad_list = []
  tad_file = open(tad_file_name, "rU")
  for line in tad_file:
    ll = line.strip().split("\t")
    if(ll[1] == "NA"): continue
    if(ll[2] == "NA"): continue
    chrom = ll[0]; start = int(ll[1]); end = int(ll[2])
    minCoord = int( floor( min(start, end) / resolution ) * resolution )
    maxCoord = int( floor( max(start, end) / resolution ) * resolution )
    tad_list.append([chrom, minCoord, maxCoord])
  tad_file.close()

  # Return objects
  return tad_list

def fetch_hic_dict(hic_file_name):

  # Structures
  hic_dict = dict()
  hic_file = open(hic_file_name, "r")

  # Reading hic_dict
  for line in hic_file:
    ll = line.strip().split("\t")
    key = ":".join(ll[:3]); value = float(ll[3])
    hic_dict[key] = value

  # Closing files and returning objects
  hic_file.close()
  return hic_dict

def fetch_hic_fc_dict(hic_1_dict, hic_2_dict, pseudocount = 1):

  # Structures
  hic_fc_dict = dict()

  # Calculating fold change
  hic_key_list = list(set().union(hic_1_dict.keys(), hic_2_dict.keys()))
  for k in hic_key_list:
    try: value1 = hic_1_dict[k] + pseudocount
    except Exception: value1 = 0 + pseudocount
    try: value2 = hic_2_dict[k] + pseudocount
    except Exception: value2 = 0 + pseudocount
    hic_fc_dict[k] = log2(value1 / value2)

  # Returning objects
  return hic_fc_dict

def chunks(vector, bins):
  index_vector = []
  if(bins == 0): return index_vector
  for i in range(0, len(vector)-1, bins):
    if(i+bins >= len(vector)-1): break
    index_vector.append(i+bins)
  return index_vector

def binned_gradient_smaller_vector(vector, length):

  # Initializing res_vec
  res_vec = [-1.0] * length
  res_vec[0] = vector[0]
  res_vec[-1] = vector[-1]

  # Getting chunks
  index_vector = chunks(vector, int(floor(float(length)/len(vector))))
  counter = 1
  for i in index_vector:
    res_vec[i] = vector[counter]
    counter += 1

  # Getting indexes to create gradient
  gradient_indexes = []
  for i in range(0,len(res_vec)):
    if(res_vec[i] >= 0): gradient_indexes.append(i)

  # Creating gradient
  for i in range(0,len(gradient_indexes)-1):
    i1 = gradient_indexes[i]
    i2 = gradient_indexes[i+1]
    v1 = res_vec[i1]
    v2 = res_vec[i2]
    if(i2-i1-1 <= 0): continue
    incr = float(v2-v1)/(i2-i1-1)
    count = 1
    for i in range(i1+1,i2):
      res_vec[i] = round(v1 + (count * incr), 4)
      count += 1

  # Returning objects
  return res_vec

def chunks2(vector, bins):
  index_vector = []
  if(bins == 0): return index_vector
  for i in range(0, len(vector)+bins, bins):
    index_vector.append(min(max(i,0),len(vector)))
  return list(OrderedDict.fromkeys(index_vector))

def binned_gradient_bigger_vector(vector, length):

  # Initializing res_vec
  res_vec = [0.0] * length

  # Calculating index of bins
  index_vector = chunks2(vector, int(floor(float(len(vector))/length)))
  if(len(index_vector)-1 != len(res_vec)):
    index_vector = index_vector[:len(res_vec)] + [index_vector[-1]]

  # Binning vector
  for i in range(0, len(index_vector)-1):
    res_vec[i] = round(float(mean(vector[index_vector[i]:index_vector[i+1]])),4)

  # Returning objects
  return res_vec

def fetch_binarized_vector(number_of_bins, signal_vector):
  
  # Creating output vector
  vector = [0.0] * number_of_bins

  # Calculating index of bins
  index_vector = chunks(signal_vector, len(signal_vector)/number_of_bins)
  if(len(index_vector)-1 != len(vector)):
    index_vector = index_vector[:len(vector)] + [index_vector[-1]]

  # Binning vector
  for i in range(0, len(index_vector)-1):
    vector[i] = round(float(mean(signal_vector[index_vector[i]:index_vector[i+1]])),4)

  # Returning objetcs
  return vector

def calculate_binned_gradient(vector, length):

  # Tads < 2 bins or number of bins < 2 not allowed
  if((len(vector) <= 1) or (length <= 1)): return None

  # Creating gradient vector
  if(len(vector) >= length): gradient_vec = binned_gradient_bigger_vector(vector, length)
  else: gradient_vec = binned_gradient_smaller_vector(vector, length)

  # Returning objects
  return gradient_vec

def update_meta_matrix(resolution, number_of_bins, meta_matrix, region1, region2, hic_1_dict, hic_2_dict):

  # Verifying if new TAD
  if(region1[1] == 0 or region1[1] == 0 or region2[1] == 0 or region2[1] == 0): return 0

  # Fetching matrix
  matrix = [[0.0] * ((region2[2]-region2[1])/resolution) for e in range(0, ((region1[2]-region1[1])/resolution))]
  counterI = 0
  chrom = region1[0]
  for i in range(region1[1], region1[2], resolution):
    counterJ = 0
    for j in range(region2[1], region2[2], resolution):
      try:
        if(i <= j): matrix[counterI][counterJ] += float(hic_2_dict[":".join([chrom,str(i),str(j)])])
        else: matrix[counterI][counterJ] += float(hic_1_dict[":".join([chrom,str(j),str(i)])])
      except Exception: pass
      counterJ += 1
    counterI += 1

  # Creating curr_meta_matrix
  curr_meta_matrix = [[0.0] * len(meta_matrix) for e in range(0,len(meta_matrix))]

  # Gradient first col
  vec = [matrix[e][0] for e in range(0,len(matrix))]
  if(vec[0] == 0): return 0
  firstCol = calculate_binned_gradient(vec, len(curr_meta_matrix))
  for j in range(0,len(curr_meta_matrix)): curr_meta_matrix[j][0] = firstCol[j]

  # Gradient last col
  vec = [matrix[e][-1] for e in range(0,len(matrix))]
  lastCol = calculate_binned_gradient(vec, len(curr_meta_matrix))
  for j in range(0,len(curr_meta_matrix)): curr_meta_matrix[j][-1] = lastCol[j]

  # Gradient by row
  for i in range(0,len(curr_meta_matrix)):
    try: curr_meta_matrix[i] = calculate_binned_gradient(matrix[i], len(curr_meta_matrix))
    except Exception: curr_meta_matrix[i] = calculate_binned_gradient([curr_meta_matrix[i][0], curr_meta_matrix[i][-1]], len(curr_meta_matrix))
    
  # Updating meta matrix
  for i in range(0,len(curr_meta_matrix)):
    for j in range(0,len(curr_meta_matrix[i])): meta_matrix[i][j] += curr_meta_matrix[i][j]

def fetch_binarized_matrix(resolution, number_of_bins, tad_1_list, tad_2_list, hic_1_dict, hic_2_dict):
  
  # Creating structures
  meta_matrix = [[0.0] * number_of_bins for e in range(0, number_of_bins)]
  
  # Iterating on regions
  for i in range(0,len(tad_1_list)):

      # Regions
      region1 = tad_1_list[i]
      region2 = tad_2_list[i]

      # Updating matrix
      update_meta_matrix(resolution, number_of_bins, meta_matrix, region1, region2, hic_1_dict, hic_2_dict)

  # Returning objetcs
  return meta_matrix

def write_meta_matrix(meta_matrix, output_file_name):

  # Writing matrix
  output_file = open(output_file_name, "w")
  for i in range(0,len(meta_matrix)):
    vec = []
    for j in range(0,len(meta_matrix)):
      if(i == j): vec.append("NA")
      else: vec.append(str(meta_matrix[i][j]))
    output_file.write("\t".join(vec)+"\n")
  output_file.close()

def create_meta_tad_table(resolution, number_of_bins, tad_1_file_name, tad_2_file_name, hic_1_file_name, hic_2_file_name, output_file_name):

  # Fetch TAD lists
  tadList1 = fetch_tad_list(resolution, tad_1_file_name)
  tadList2 = fetch_tad_list(resolution, tad_2_file_name)

  # Fetch hic_dict
  hic_1_dict = fetch_hic_dict(hic_1_file_name)
  hic_2_dict = fetch_hic_dict(hic_2_file_name)

  # Creating meta matrix
  meta_matrix = fetch_binarized_matrix(resolution, number_of_bins, tadList1, tadList2, hic_1_dict, hic_2_dict)

  # Writing meta matrix
  write_meta_matrix(meta_matrix, output_file_name)

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
create_meta_tad_table(resolution, number_of_bins, tad_1_file_name, tad_2_file_name, hic_1_file_name, hic_2_file_name, output_file_name)


